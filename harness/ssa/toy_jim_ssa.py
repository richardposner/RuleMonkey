#!/usr/bin/env python3
"""
Custom Gillespie SSA simulator for the toy_jim model.

Ground-truth simulator for verifying RM and NFsim correctness on this model.
Tracks all molecules, bonds, complexes, and phosphorylation states explicitly.

Model: L(r), R(l,r,a), A(r,k), K(a,Y~U/P)
  1000 each, t_end=100, n_steps=100

Rules (BNGL):
  R1:  L(r) + R(l,r)       <-> L(r!1).R(l!1,r)                   kpL, kmL
  R2:  L(r!1).R(l!1,r) + L(r!1).R(l!1,r) <-> dimer               kpD, kmD
  R3:  A(r) + R(a)          <-> A(r!1).R(a!1)                     kpA, kmA
  R4:  A(k) + K(a)          <-> A(k!1).K(a!1)                     kpK, kmK
  R5:  K(Y~U).K(Y~U)        ->  K(Y~U).K(Y~P)                    pK
  R6:  K(Y~P).K(Y~U)        ->  K(Y~P).K(Y~P)                    pKs
  R7:  R(a!1).A(r!1,k!2).K(a!2,Y~P) -> R(a!1).A(r!1,k!2).K(a!2,Y~U)  dM
  R8:  K(a,Y~P)             ->  K(a,Y~U)                          dC

Usage:
  python3 toy_jim_ssa.py                         # 100 reps, write ensemble
  python3 toy_jim_ssa.py --reps 5 --seed 42      # 5 reps, fixed seed
  python3 toy_jim_ssa.py --gdat                   # single rep, print gdat to stdout
"""

import sys, os, math, random, argparse, time as timepkg
from collections import defaultdict

# ── Model parameters ──────────────────────────────────────────────────────────

N_L = 1000
N_R = 1000
N_A = 1000
N_K = 1000
N_MOLS = N_L + N_R + N_A + N_K  # 4000

kpL = 0.0001
kmL = 0.1
kpD = 0.001    # already /2 from toy.in for BNG symmetry
kmD = 0.1
kpA = 0.0001
kmA = 0.1
kpK = 0.0001
kmK = 0.1
pK  = 1.0
pKs = 0.01
dM  = 1.0
dC  = 10.0

T_END   = 100.0
N_STEPS = 100

# Molecule index ranges
L_START = 0
L_END   = N_L
R_START = N_L
R_END   = N_L + N_R
A_START = N_L + N_R
A_END   = N_L + N_R + N_A
K_START = N_L + N_R + N_A
K_END   = N_MOLS

# Sites per molecule type
# L: site 0 = r
# R: site 0 = l, site 1 = r, site 2 = a
# A: site 0 = r, site 1 = k
# K: site 0 = a  (Y is internal state, not a bond site)

OBS_NAMES = ["RecDim_1", "Rec_A_1", "Rec_K_1", "Rec_Kp_1", "RecDim_Kp_1",
             "L_tot_1", "A_tot_1", "K_tot_1", "R_tot_1"]


class Simulator:
    __slots__ = (
        "bp",          # bond partners: bp[mol][site] = (other_mol, other_site) or None
        "k_phos",      # k_phos[k_mol_idx] = True if phosphorylated
        "cx_of",       # mol -> complex_id
        "cx_mols",     # cx_id -> set of mol_ids
        "next_cx",
        "t",
        # Aggregate counts maintained incrementally
        "n_free_L",       # L with r free (= free L molecules, since monovalent)
        "n_R_lr_free",    # R with both l and r free
        "n_LR_mono",      # L-R bonds where R.r is free (eligible for dimerization)
        "n_RR_bonds",     # number of R-R dimer bonds
        "n_free_Ar",      # A with r free
        "n_free_Ra",      # R with a free
        "n_free_Ak",      # A with k free
        "n_free_Ka",      # K with a free
        "n_free_Ka_YP",   # K with a free and Y~P
        "n_RAK_YP",       # R-A-K(Y~P) chains (count of such chains)
        # Per-complex K state tracking for transphosphorylation
        "cx_nKU",         # cx_id -> count of K(Y~U) in complex
        "cx_nKP",         # cx_id -> count of K(Y~P) in complex
        "sum_KU_pairs",   # Σ_cx n_KU*(n_KU-1)/2  (rule 5 propensity factor)
        "sum_KP_KU",      # Σ_cx n_KP*n_KU          (rule 6 propensity factor)
        # Bond lists for uniform selection
        "lr_bonds",       # list of (L_mol, R_mol) where L.r-R.l bond exists and R.r free
        "rr_bonds",       # list of (R1_mol, R2_mol) where R1.r-R2.r bond exists
        "ar_bonds",       # list of (A_mol, R_mol) where A.r-R.a bond exists
        "ak_bonds",       # list of (A_mol, K_mol) where A.k-K.a bond exists
    )

    def __init__(self, seed=None):
        if seed is not None:
            random.seed(seed)

        # Bond partners: None = free
        self.bp = []
        for i in range(N_MOLS):
            if i < L_END:
                self.bp.append([None])          # L: 1 site
            elif i < R_END:
                self.bp.append([None, None, None])  # R: 3 sites
            elif i < A_END:
                self.bp.append([None, None])    # A: 2 sites
            else:
                self.bp.append([None])          # K: 1 site (a)

        self.k_phos = [False] * N_K  # index within K range: k_phos[k - K_START]

        # Each molecule is its own complex initially
        self.cx_of = list(range(N_MOLS))
        self.cx_mols = {i: {i} for i in range(N_MOLS)}
        self.next_cx = N_MOLS

        self.cx_nKU = defaultdict(int)
        self.cx_nKP = defaultdict(int)
        for k in range(K_START, K_END):
            self.cx_nKU[k] = 1  # all K start as Y~U

        # Aggregate counts
        self.n_free_L = N_L
        self.n_R_lr_free = N_R    # all R start with l, r, a free
        self.n_LR_mono = 0
        self.n_RR_bonds = 0
        self.n_free_Ar = N_A
        self.n_free_Ra = N_R
        self.n_free_Ak = N_A
        self.n_free_Ka = N_K
        self.n_free_Ka_YP = 0
        self.n_RAK_YP = 0

        self.sum_KU_pairs = 0
        self.sum_KP_KU = 0

        self.lr_bonds = []
        self.rr_bonds = []
        self.ar_bonds = []
        self.ak_bonds = []

        self.t = 0.0

    # ── helpers ───────────────────────────────────────────────────────────

    def _mol_type(self, m):
        if m < L_END: return 'L'
        if m < R_END: return 'R'
        if m < A_END: return 'A'
        return 'K'

    def _is_K_phos(self, k):
        return self.k_phos[k - K_START]

    def _set_K_phos(self, k, val):
        self.k_phos[k - K_START] = val

    def _merge_cx(self, m1, m2):
        """Merge complexes of m1 and m2. Returns the surviving cx_id."""
        c1, c2 = self.cx_of[m1], self.cx_of[m2]
        if c1 == c2:
            return c1
        # merge smaller into larger
        if len(self.cx_mols[c1]) < len(self.cx_mols[c2]):
            c1, c2 = c2, c1
        # Update sum_KU_pairs and sum_KP_KU: remove old contributions, add merged
        nku1, nkp1 = self.cx_nKU[c1], self.cx_nKP[c1]
        nku2, nkp2 = self.cx_nKU[c2], self.cx_nKP[c2]
        self.sum_KU_pairs -= nku1 * (nku1 - 1) // 2 + nku2 * (nku2 - 1) // 2
        self.sum_KP_KU -= nkp1 * nku1 + nkp2 * nku2
        nku = nku1 + nku2
        nkp = nkp1 + nkp2
        self.sum_KU_pairs += nku * (nku - 1) // 2
        self.sum_KP_KU += nkp * nku

        for m in self.cx_mols[c2]:
            self.cx_of[m] = c1
        self.cx_mols[c1] |= self.cx_mols[c2]
        self.cx_nKU[c1] = nku
        self.cx_nKP[c1] = nkp
        del self.cx_mols[c2]
        del self.cx_nKU[c2]
        del self.cx_nKP[c2]
        return c1

    def _split_check(self, m1, m2):
        """After removing a bond between m1 and m2, check if they're still connected.
        If not, split the complex."""
        cx = self.cx_of[m1]
        # BFS from m1, see if we reach m2
        visited = {m1}
        stack = [m1]
        while stack:
            cur = stack.pop()
            for partner in self.bp[cur]:
                if partner is None:
                    continue
                nb = partner[0]
                if nb == m2:
                    return  # still connected
                if nb not in visited:
                    visited.add(nb)
                    stack.append(nb)
        # m2 not reachable from m1 → split
        # visited = component of m1
        # Everything else in cx_mols[cx] that's not in visited = component of m2
        comp2 = self.cx_mols[cx] - visited

        # Remove old contributions
        nku_old, nkp_old = self.cx_nKU[cx], self.cx_nKP[cx]
        self.sum_KU_pairs -= nku_old * (nku_old - 1) // 2
        self.sum_KP_KU -= nkp_old * nku_old

        # Compute K counts for component 2
        nku2 = sum(1 for m in comp2 if m >= K_START and not self._is_K_phos(m))
        nkp2 = sum(1 for m in comp2 if m >= K_START and self._is_K_phos(m))
        nku1 = nku_old - nku2
        nkp1 = nkp_old - nkp2

        # Create new complex for comp2
        ncx = self.next_cx
        self.next_cx += 1
        for m in comp2:
            self.cx_of[m] = ncx
        self.cx_mols[cx] -= comp2
        self.cx_mols[ncx] = comp2

        self.cx_nKU[cx] = nku1
        self.cx_nKP[cx] = nkp1
        self.cx_nKU[ncx] = nku2
        self.cx_nKP[ncx] = nkp2

        # Add new contributions
        self.sum_KU_pairs += nku1 * (nku1 - 1) // 2 + nku2 * (nku2 - 1) // 2
        self.sum_KP_KU += nkp1 * nku1 + nkp2 * nku2

    def _remove_from_list(self, lst, item):
        """Remove item from list by swapping with last element."""
        try:
            idx = lst.index(item)
        except ValueError:
            # Try reversed tuple
            idx = lst.index((item[1], item[0]))
        lst[idx] = lst[-1]
        lst.pop()

    # ── propensities ─────────────────────────────────────────────────────

    def propensities(self):
        # R1 fwd: L(r) + R(l,r) -> L(r!1).R(l!1,r)
        # Free L is always singleton, so inter-complex automatic
        a1f = kpL * self.n_free_L * self.n_R_lr_free

        # R1 rev: L(r!1).R(l!1,r) -> L(r) + R(l,r)
        # L-R bond where R.r is free
        a1r = kmL * self.n_LR_mono

        # R2 fwd: dimerization of LR monomers (identical reactant patterns)
        # Two LR monomers can never be in the same complex (see analysis in code comments)
        # BNG convention: propensity = kpD * C(n, 2) = kpD * n*(n-1)/2
        n = self.n_LR_mono
        a2f = kpD * n * (n - 1) / 2

        # R2 rev: dimer dissociation (R.r-R.r bond break)
        a2r = kmD * self.n_RR_bonds

        # R3 fwd: A(r) + R(a) -> A(r!1).R(a!1)
        # A with r free is either singleton or A-K pair; R with a free is in some complex
        # They cannot be in the same complex (A.r is the only link from A-side to R-side)
        a3f = kpA * self.n_free_Ar * self.n_free_Ra

        # R3 rev: A(r!1).R(a!1) -> A(r) + R(a)
        a3r = kmA * len(self.ar_bonds)

        # R4 fwd: A(k) + K(a) -> A(k!1).K(a!1)
        # K with a free is always singleton; A with k free may be in any complex
        a4f = kpK * self.n_free_Ak * self.n_free_Ka

        # R4 rev: A(k!1).K(a!1) -> A(k) + K(a)
        a4r = kmK * len(self.ak_bonds)

        # R5: K(Y~U).K(Y~U) -> K(Y~U).K(Y~P)
        # Disjoint same-complex pattern; symmetry_factor=1 in BNG
        # propensity = pK * Σ_cx n_KU*(n_KU-1) [both ordered embeddings count, sf=1]
        # Wait: sf=1 means propensity = a_total * rate * 1
        # a_total = Σ_cx n_KU*(n_KU-1) [seeding from each K, finding (n_KU-1) partners]
        # = 2 * Σ_cx C(n_KU, 2)
        # So propensity = 2 * sum_KU_pairs * pK
        a5 = pK * 2 * self.sum_KU_pairs

        # R6: K(Y~P).K(Y~U) -> K(Y~P).K(Y~P)
        # Disjoint same-complex pattern; no symmetry (different patterns)
        # propensity = pKs * Σ_cx n_KP * n_KU
        a6 = pKs * self.sum_KP_KU

        # R7: R(a!1).A(r!1,k!2).K(a!2,Y~P) -> R(a!1).A(r!1,k!2).K(a!2,Y~U)
        a7 = dM * self.n_RAK_YP

        # R8: K(a,Y~P) -> K(a,Y~U)
        a8 = dC * self.n_free_Ka_YP

        return (a1f, a1r, a2f, a2r, a3f, a3r, a4f, a4r, a5, a6, a7, a8)

    # ── reaction firing ──────────────────────────────────────────────────

    def fire(self, channel):
        if channel == 0:
            self._fire_R1_fwd()
        elif channel == 1:
            self._fire_R1_rev()
        elif channel == 2:
            self._fire_R2_fwd()
        elif channel == 3:
            self._fire_R2_rev()
        elif channel == 4:
            self._fire_R3_fwd()
        elif channel == 5:
            self._fire_R3_rev()
        elif channel == 6:
            self._fire_R4_fwd()
        elif channel == 7:
            self._fire_R4_rev()
        elif channel == 8:
            self._fire_R5()
        elif channel == 9:
            self._fire_R6()
        elif channel == 10:
            self._fire_R7()
        elif channel == 11:
            self._fire_R8()

    def _pick_free_L(self):
        """Pick a random free L molecule."""
        candidates = [m for m in range(L_START, L_END) if self.bp[m][0] is None]
        return random.choice(candidates)

    def _pick_R_lr_free(self):
        """Pick a random R with both l and r free."""
        candidates = [m for m in range(R_START, R_END)
                      if self.bp[m][0] is None and self.bp[m][1] is None]
        return random.choice(candidates)

    def _fire_R1_fwd(self):
        """L(r) + R(l,r) -> L(r!1).R(l!1,r)"""
        L = self._pick_free_L()
        R = self._pick_R_lr_free()
        # Create bond L.r(0) - R.l(0)
        self.bp[L][0] = (R, 0)
        self.bp[R][0] = (L, 0)
        # Update counts
        self.n_free_L -= 1
        self.n_R_lr_free -= 1
        self.n_LR_mono += 1
        self.lr_bonds.append((L, R))
        # R.a status unchanged, so n_free_Ra unchanged
        # Merge complexes (L was singleton)
        self._merge_cx(L, R)

    def _fire_R1_rev(self):
        """Unbind L-R bond where R.r is free."""
        idx = random.randrange(len(self.lr_bonds))
        L, R = self.lr_bonds[idx]
        # Verify R.r is still free (should be by definition of lr_bonds)
        assert self.bp[R][1] is None, "lr_bonds entry has R.r bound!"
        # Remove bond
        self.bp[L][0] = None
        self.bp[R][0] = None
        # Update bond list
        self.lr_bonds[idx] = self.lr_bonds[-1]
        self.lr_bonds.pop()
        # Update counts
        self.n_free_L += 1
        self.n_R_lr_free += 1
        self.n_LR_mono -= 1
        # Split complex
        self._split_check(L, R)

    def _fire_R2_fwd(self):
        """Dimerize two LR monomers."""
        # Pick 2 distinct lr_bonds entries
        n = len(self.lr_bonds)
        i = random.randrange(n)
        j = random.randrange(n - 1)
        if j >= i:
            j += 1
        L1, R1 = self.lr_bonds[i]
        L2, R2 = self.lr_bonds[j]
        # Create R1.r(1) - R2.r(1) bond
        self.bp[R1][1] = (R2, 1)
        self.bp[R2][1] = (R1, 1)
        # Remove both from lr_bonds (they're no longer monomeric)
        # Remove higher index first to avoid shifting
        if i > j:
            self.lr_bonds[i] = self.lr_bonds[-1]
            self.lr_bonds.pop()
            self.lr_bonds[j] = self.lr_bonds[-1]
            self.lr_bonds.pop()
        else:
            self.lr_bonds[j] = self.lr_bonds[-1]
            self.lr_bonds.pop()
            self.lr_bonds[i] = self.lr_bonds[-1]
            self.lr_bonds.pop()
        # Update counts
        self.n_LR_mono -= 2
        self.n_RR_bonds += 1
        self.rr_bonds.append((R1, R2))
        # Merge complexes
        self._merge_cx(R1, R2)

    def _fire_R2_rev(self):
        """Break R-R dimer bond."""
        idx = random.randrange(len(self.rr_bonds))
        R1, R2 = self.rr_bonds[idx]
        # Remove bond
        self.bp[R1][1] = None
        self.bp[R2][1] = None
        self.rr_bonds[idx] = self.rr_bonds[-1]
        self.rr_bonds.pop()
        self.n_RR_bonds -= 1
        # Both R's now have r free. If they still have L bound (l site), they become LR monomers
        if self.bp[R1][0] is not None:
            self.n_LR_mono += 1
            L1 = self.bp[R1][0][0]
            self.lr_bonds.append((L1, R1))
        if self.bp[R2][0] is not None:
            self.n_LR_mono += 1
            L2 = self.bp[R2][0][0]
            self.lr_bonds.append((L2, R2))
        # Split complex
        self._split_check(R1, R2)

    def _fire_R3_fwd(self):
        """A(r) + R(a) -> A(r!1).R(a!1)"""
        # Pick random A with r free
        candidates_A = [m for m in range(A_START, A_END) if self.bp[m][0] is None]
        A = random.choice(candidates_A)
        # Pick random R with a free
        candidates_R = [m for m in range(R_START, R_END) if self.bp[m][2] is None]
        R = random.choice(candidates_R)
        # Create bond A.r(0) - R.a(2)
        self.bp[A][0] = (R, 2)
        self.bp[R][2] = (A, 0)
        self.n_free_Ar -= 1
        self.n_free_Ra -= 1
        self.ar_bonds.append((A, R))
        # Check if this creates an R-A-K(Y~P) chain
        if self.bp[A][1] is not None:
            K = self.bp[A][1][0]
            if self._is_K_phos(K):
                self.n_RAK_YP += 1
        # Merge complexes
        self._merge_cx(A, R)

    def _fire_R3_rev(self):
        """Break A.r-R.a bond."""
        idx = random.randrange(len(self.ar_bonds))
        A, R = self.ar_bonds[idx]
        # Check if we're breaking an R-A-K(Y~P) chain
        if self.bp[A][1] is not None:
            K = self.bp[A][1][0]
            if self._is_K_phos(K):
                self.n_RAK_YP -= 1
        # Remove bond
        self.bp[A][0] = None
        self.bp[R][2] = None
        self.ar_bonds[idx] = self.ar_bonds[-1]
        self.ar_bonds.pop()
        self.n_free_Ar += 1
        self.n_free_Ra += 1
        # Split complex
        self._split_check(A, R)

    def _fire_R4_fwd(self):
        """A(k) + K(a) -> A(k!1).K(a!1)"""
        candidates_A = [m for m in range(A_START, A_END) if self.bp[m][1] is None]
        A = random.choice(candidates_A)
        candidates_K = [m for m in range(K_START, K_END) if self.bp[m][0] is None]
        K = random.choice(candidates_K)
        # Create bond A.k(1) - K.a(0)
        self.bp[A][1] = (K, 0)
        self.bp[K][0] = (A, 1)
        self.n_free_Ak -= 1
        self.n_free_Ka -= 1
        if self._is_K_phos(K):
            self.n_free_Ka_YP -= 1
        self.ak_bonds.append((A, K))
        # Check if this creates an R-A-K chain (A must be bound to R via A.r)
        if self.bp[A][0] is not None:
            if self._is_K_phos(K):
                self.n_RAK_YP += 1
        # Merge complexes
        self._merge_cx(A, K)

    def _fire_R4_rev(self):
        """Break A.k-K.a bond."""
        idx = random.randrange(len(self.ak_bonds))
        A, K = self.ak_bonds[idx]
        # Check if we're breaking an R-A-K(Y~P) chain
        if self.bp[A][0] is not None and self._is_K_phos(K):
            self.n_RAK_YP -= 1
        # Remove bond
        self.bp[A][1] = None
        self.bp[K][0] = None
        self.ak_bonds[idx] = self.ak_bonds[-1]
        self.ak_bonds.pop()
        self.n_free_Ak += 1
        self.n_free_Ka += 1
        if self._is_K_phos(K):
            self.n_free_Ka_YP += 1
        # Split complex
        self._split_check(A, K)

    def _fire_R5(self):
        """K(Y~U).K(Y~U) -> K(Y~U).K(Y~P) — inactive transphosphorylation.
        symmetry_factor=1: each ordered embedding is a distinct reaction.
        We pick uniformly among all ordered (seed, partner) pairs."""
        # Total ordered pairs = 2 * sum_KU_pairs = Σ n*(n-1)
        # Pick a complex weighted by n*(n-1), then pick an ordered pair
        target = random.random() * 2 * self.sum_KU_pairs
        cumul = 0.0
        for cx_id, mols in self.cx_mols.items():
            nku = self.cx_nKU.get(cx_id, 0)
            w = nku * (nku - 1)  # ordered pairs
            if w <= 0:
                continue
            cumul += w
            if target < cumul:
                # Pick an ordered pair of K(Y~U) in this complex
                kus = [m for m in mols if m >= K_START and not self._is_K_phos(m)]
                # Pick the K to phosphorylate (the "second" one)
                victim = random.choice(kus)
                self._phosphorylate(victim)
                return
        # Fallback (shouldn't happen if sum_KU_pairs > 0)

    def _fire_R6(self):
        """K(Y~P).K(Y~U) -> K(Y~P).K(Y~P) — active transphosphorylation."""
        # Each (KP, KU) ordered pair is a distinct reaction (no symmetry)
        target = random.random() * self.sum_KP_KU
        cumul = 0.0
        for cx_id, mols in self.cx_mols.items():
            nkp = self.cx_nKP.get(cx_id, 0)
            nku = self.cx_nKU.get(cx_id, 0)
            w = nkp * nku
            if w <= 0:
                continue
            cumul += w
            if target < cumul:
                # Pick a random K(Y~U) in this complex to phosphorylate
                kus = [m for m in mols if m >= K_START and not self._is_K_phos(m)]
                victim = random.choice(kus)
                self._phosphorylate(victim)
                return

    def _phosphorylate(self, k):
        """Change K from Y~U to Y~P."""
        cx = self.cx_of[k]
        nku = self.cx_nKU[cx]
        nkp = self.cx_nKP[cx]
        # Remove old contributions
        self.sum_KU_pairs -= nku * (nku - 1) // 2
        self.sum_KP_KU -= nkp * nku
        # Update state
        self._set_K_phos(k, True)
        self.cx_nKU[cx] -= 1
        self.cx_nKP[cx] += 1
        nku -= 1
        nkp += 1
        # Add new contributions
        self.sum_KU_pairs += nku * (nku - 1) // 2
        self.sum_KP_KU += nkp * nku
        # Update observable-related counts
        if self.bp[k][0] is None:
            # K.a is free → update free_Ka_YP
            self.n_free_Ka_YP += 1
        else:
            # K is bound to A via K.a-A.k → check if R-A-K chain exists
            A = self.bp[k][0][0]
            if self.bp[A][0] is not None:
                # A.r bound to R → this is now an R-A-K(Y~P) chain
                self.n_RAK_YP += 1

    def _dephosphorylate(self, k):
        """Change K from Y~P to Y~U."""
        cx = self.cx_of[k]
        nku = self.cx_nKU[cx]
        nkp = self.cx_nKP[cx]
        # Remove old contributions
        self.sum_KU_pairs -= nku * (nku - 1) // 2
        self.sum_KP_KU -= nkp * nku
        # Update state
        self._set_K_phos(k, False)
        self.cx_nKU[cx] += 1
        self.cx_nKP[cx] -= 1
        nku += 1
        nkp -= 1
        # Add new contributions
        self.sum_KU_pairs += nku * (nku - 1) // 2
        self.sum_KP_KU += nkp * nku
        # Update observable-related counts
        if self.bp[k][0] is None:
            self.n_free_Ka_YP -= 1
        else:
            A = self.bp[k][0][0]
            if self.bp[A][0] is not None:
                self.n_RAK_YP -= 1

    def _fire_R7(self):
        """R(a!1).A(r!1,k!2).K(a!2,Y~P) -> dephosphorylate membrane K."""
        # Pick a random R-A-K(Y~P) chain
        # Enumerate all such chains
        chains = []
        for A, R in self.ar_bonds:
            if self.bp[A][1] is not None:
                K = self.bp[A][1][0]
                if self._is_K_phos(K):
                    chains.append(K)
        K = random.choice(chains)
        self._dephosphorylate(K)

    def _fire_R8(self):
        """K(a,Y~P) -> K(a,Y~U) — cytosol dephosphorylation."""
        # K with a free and Y~P
        candidates = [m for m in range(K_START, K_END)
                      if self.bp[m][0] is None and self._is_K_phos(m)]
        K = random.choice(candidates)
        self._dephosphorylate(K)

    # ── observables ──────────────────────────────────────────────────────

    def observe(self):
        # RecDim_1: Molecules R(r!+) — count of R with r bound
        rec_dim = sum(1 for m in range(R_START, R_END) if self.bp[m][1] is not None)

        # Rec_A_1: Molecules R(a!1).A(r!1) — count of R-A bonds
        rec_a = len(self.ar_bonds)

        # Rec_K_1: Molecules R(a!1).A(r!1,k!2).K(a!2)
        # = R-A-K chains (any K phosphorylation state)
        rec_k = 0
        for A, R in self.ar_bonds:
            if self.bp[A][1] is not None:
                rec_k += 1

        # Rec_Kp_1: Molecules R(a!1).A(r!1,k!2).K(a!2,Y~P)
        rec_kp = self.n_RAK_YP

        # RecDim_Kp_1: Molecules R().R(a!1).A(r!1,k!2).K(a!2,Y~P)
        # For each R-A-K(Y~P) chain, count (n_R_in_complex - 1)
        rec_dim_kp = 0
        for A, R in self.ar_bonds:
            if self.bp[A][1] is not None:
                K = self.bp[A][1][0]
                if self._is_K_phos(K):
                    cx = self.cx_of[R]
                    n_R_in_cx = sum(1 for m in self.cx_mols[cx] if R_START <= m < R_END)
                    rec_dim_kp += (n_R_in_cx - 1)

        return [rec_dim, rec_a, rec_k, rec_kp, rec_dim_kp,
                N_L, N_A, N_K, N_R]

    # ── main loop ────────────────────────────────────────────────────────

    def run(self):
        dt = T_END / N_STEPS
        results = [self.observe()]  # t=0

        next_out = dt
        out_idx = 1
        event_count = 0

        while out_idx <= N_STEPS:
            ps = self.propensities()
            total = sum(ps)
            if total <= 0:
                while out_idx <= N_STEPS:
                    results.append(self.observe())
                    out_idx += 1
                break

            tau = -math.log(random.random()) / total

            while out_idx <= N_STEPS and self.t + tau >= next_out:
                results.append(self.observe())
                next_out += dt
                out_idx += 1

            self.t += tau
            if self.t >= T_END:
                break

            # Select channel
            r = random.random() * total
            cumul = 0.0
            channel = 0
            for i, p in enumerate(ps):
                cumul += p
                if r < cumul:
                    channel = i
                    break

            self.fire(channel)
            event_count += 1

        while len(results) <= N_STEPS:
            results.append(self.observe())

        return results, event_count


# ── Ensemble generation ──────────────────────────────────────────────────────

def run_ensemble(n_reps, base_seed=None):
    import numpy as np
    n_obs = len(OBS_NAMES)
    all_data = np.zeros((n_reps, N_STEPS + 1, n_obs))

    for rep in range(n_reps):
        seed = (base_seed + rep) if base_seed is not None else None
        t0 = timepkg.time()
        sim = Simulator(seed=seed)
        rows, events = sim.run()
        elapsed = timepkg.time() - t0
        print(f"  rep {rep+1}/{n_reps}  events={events}  {elapsed:.1f}s",
              file=sys.stderr)
        for ti, row in enumerate(rows):
            all_data[rep, ti, :] = row

    return all_data


def write_ensemble(all_data, out_dir):
    import numpy as np
    os.makedirs(out_dir, exist_ok=True)

    n_reps = all_data.shape[0]
    dt = T_END / N_STEPS
    times = [dt * i for i in range(N_STEPS + 1)]

    mean = np.mean(all_data, axis=0)
    std = np.std(all_data, axis=0, ddof=1)

    with open(os.path.join(out_dir, "toy_jim.mean.tsv"), "w") as f:
        f.write("time\t" + "\t".join(OBS_NAMES) + "\n")
        for ti in range(N_STEPS + 1):
            vals = "\t".join(f"{mean[ti, oi]:.6g}" for oi in range(len(OBS_NAMES)))
            f.write(f"{times[ti]:.6g}\t{vals}\n")

    with open(os.path.join(out_dir, "toy_jim.std.tsv"), "w") as f:
        f.write("time\t" + "\t".join(OBS_NAMES) + "\n")
        for ti in range(N_STEPS + 1):
            vals = "\t".join(f"{std[ti, oi]:.6g}" for oi in range(len(OBS_NAMES)))
            f.write(f"{times[ti]:.6g}\t{vals}\n")

    with open(os.path.join(out_dir, "toy_jim.tint.tsv"), "w") as f:
        f.write("observable\tI_mean\tI_std\tn_reps\n")
        for oi, name in enumerate(OBS_NAMES):
            integrals = np.zeros(n_reps)
            for r in range(n_reps):
                for ti in range(N_STEPS):
                    integrals[r] += 0.5 * (all_data[r, ti, oi] + all_data[r, ti+1, oi]) * dt
            f.write(f"{name}\t{np.mean(integrals):.6g}\t{np.std(integrals, ddof=1):.6g}"
                    f"\t{n_reps}\n")

    rep_dir = os.path.join(
        os.path.dirname(out_dir), "replicates", "toy_jim")
    os.makedirs(rep_dir, exist_ok=True)
    for r in range(n_reps):
        tag = f"rep_{r:03d}"
        gdat = os.path.join(rep_dir, f"{tag}_toy_jim.gdat")
        with open(gdat, "w") as f:
            hdr = "#" + " ".join(f"{c:>15s}" for c in ["time"] + list(OBS_NAMES))
            f.write(hdr + "\n")
            for ti in range(N_STEPS + 1):
                vals = [f"{times[ti]:.8e}"] + [
                    f"{all_data[r, ti, oi]:.8e}" for oi in range(len(OBS_NAMES))]
                f.write("\t".join(vals) + "\n")
        open(os.path.join(rep_dir, f"{tag}.done"), "w").close()

    print(f"Wrote ensemble to {out_dir}/", file=sys.stderr)
    print(f"Wrote {n_reps} per-rep gdat files to {rep_dir}/", file=sys.stderr)


def print_gdat(rows):
    dt = T_END / N_STEPS
    header = "# time\t" + "\t".join(OBS_NAMES)
    print(header)
    for ti, row in enumerate(rows):
        t = dt * ti
        vals = "\t".join(str(v) for v in row)
        print(f"{t}\t{vals}")


def main():
    ap = argparse.ArgumentParser(description="toy_jim Gillespie SSA")
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--output", type=str, default=None,
                    help="Output directory for ensemble files")
    ap.add_argument("--gdat", action="store_true",
                    help="Single rep, print gdat to stdout")
    args = ap.parse_args()

    if args.gdat:
        sim = Simulator(seed=args.seed)
        rows, events = sim.run()
        print_gdat(rows)
        print(f"Events: {events}", file=sys.stderr)
        return

    out_dir = args.output or os.path.join(
        os.path.dirname(__file__), "..", "models", "nfsim_reference", "ensemble")

    print(f"Running {args.reps} reps with seed={args.seed}", file=sys.stderr)
    t0 = timepkg.time()
    all_data = run_ensemble(args.reps, base_seed=args.seed)
    elapsed = timepkg.time() - t0
    print(f"Total: {elapsed:.1f}s ({elapsed/args.reps:.1f}s/rep)", file=sys.stderr)

    write_ensemble(all_data, out_dir)


if __name__ == "__main__":
    main()
