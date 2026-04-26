#!/usr/bin/env python3
"""
Custom Gillespie SSA simulator for the rm_tlbr_rings model.

This purpose-built simulator exists because NFsim has two confirmed bugs
that prevent it from generating correct reference data for this model:

1. NFsim's -bscb flag does not enforce product-side "+" for unimolecular
   reverse rules (RuleWorld/nfsim#61).  Bond-breaking rules fire even
   when the products remain in the same complex (ring bonds).

2. NFsim under-counts embeddings for disjoint (connectedTo) patterns
   when molecules have symmetric components (RuleWorld/nfsim#62).
   The ring-closure rule L(r).R(l) is a disjoint pattern; NFsim finds
   fewer valid L-R intra-complex pairs than actually exist.

Other ring models (BLBR, A+B rings) avoid both bugs because they use
fully-connected bond patterns for ring closure.

The simulator implements strict BNGL semantics:
- Bimolecular rules enforce inter-complex only (-bscb).
- Bond-breaking with "+" products fires only when products separate (bridge bonds).
- Bond-breaking with ensureConnected fires only when products stay connected (ring bonds).
- Embedding counts are exact (no -utl approximation).
- All unbinding rates are equal (km1=km2=km3=0.01), so every bond breaks at the
  same rate regardless of category.

Model: rm_tlbr_rings — L(r,r,r) + R(l,l), 4200 L, 300 R, t_end=1000, n_steps=100

Usage:
  python3 rm_tlbr_rings_ssa.py                     # 100 reps, write ensemble
  python3 rm_tlbr_rings_ssa.py --reps 5 --seed 42  # 5 reps, fixed seed
  python3 rm_tlbr_rings_ssa.py --gdat               # single rep, print gdat to stdout
"""

import sys, os, math, random, argparse, time as timepkg
from collections import defaultdict

# ── Model parameters ──────────────────────────────────────────────────────────

LIG_TOT = 4200
REC_TOT = 300
KP1 = 3.0e-7   # first capture (L with 0 bonds + R, inter-complex)
KP2 = 3.0e-4   # crosslink    (L with 1+ bonds + R, inter-complex)
KP3 = 3.0e-4   # ring closure (L + R, intra-complex)  [== KP2]
KM  = 0.01     # all unbinding (km1 == km2 == km3)
T_END = 1000.0
N_STEPS = 100
N_MOLS = LIG_TOT + REC_TOT

OBS_NAMES = (
    ["Ligfree_1", "Ligbnd1_1", "Ligbnd2_1", "Ligbnd3_1"]
    + [f"Size_{n}" for n in range(1, REC_TOT + 1)]
)

# ── Simulator ─────────────────────────────────────────────────────────────────

class Simulator:
    __slots__ = (
        "mol_type", "nsites", "bp", "nbonds", "cx_of", "cx_mols",
        "cx_freeL", "cx_freeR", "bond_list", "total_bonds",
        "n_freeL_L0", "n_freeL_L1", "n_freeR", "intra_pairs",
        "next_cx", "t",
    )

    def __init__(self, seed=None):
        if seed is not None:
            random.seed(seed)

        # mol 0..LIG_TOT-1 = L (3 sites), LIG_TOT..N_MOLS-1 = R (2 sites)
        self.mol_type = [0] * LIG_TOT + [1] * REC_TOT   # 0=L, 1=R
        self.nsites   = [3] * LIG_TOT + [2] * REC_TOT
        # bond partner: bp[mol][site] = (other_mol, other_site) or None
        self.bp = [
            [None, None, None] if i < LIG_TOT else [None, None]
            for i in range(N_MOLS)
        ]
        self.nbonds = [0] * N_MOLS            # bonds on each molecule

        # complexes  (initially every molecule is its own complex)
        self.cx_of   = list(range(N_MOLS))     # mol -> cx_id
        self.cx_mols = {i: {i} for i in range(N_MOLS)}
        self.cx_freeL = defaultdict(int)        # cx -> free L r-sites in cx
        self.cx_freeR = defaultdict(int)
        for i in range(LIG_TOT):
            self.cx_freeL[i] = 3
        for i in range(LIG_TOT, N_MOLS):
            self.cx_freeR[i] = 2
        self.next_cx = N_MOLS

        # bond list for uniform sampling on unbinding
        self.bond_list = []   # [(molL, siteL, molR, siteR), ...]
        self.total_bonds = 0

        # aggregate counts
        self.n_freeL_L0 = 3 * LIG_TOT   # free r-sites on L with 0 bonds
        self.n_freeL_L1 = 0              # free r-sites on L with 1+ bonds
        self.n_freeR    = 2 * REC_TOT
        self.intra_pairs = 0             # Σ_cx (freeL_cx × freeR_cx)

        self.t = 0.0

    # ── bond operations ───────────────────────────────────────────────────

    def _add_bond(self, mL, sL, mR, sR):
        cxL = self.cx_of[mL]
        cxR = self.cx_of[mR]
        same = (cxL == cxR)

        # remove old intra contribution
        old_i = self.cx_freeL[cxL] * self.cx_freeR[cxL]
        if not same:
            old_i += self.cx_freeL[cxR] * self.cx_freeR[cxR]

        # update molecule-level counts
        ob = self.nbonds[mL]
        self.nbonds[mL] += 1
        self.nbonds[mR] += 1

        # update aggregate free-site counts
        free_after = (3 - self.nbonds[mL])   # free r-sites on mL after
        if ob == 0:
            self.n_freeL_L0 -= (free_after + 1)  # was free_after+1 in L0
            self.n_freeL_L1 += free_after
        else:
            self.n_freeL_L1 -= 1
        self.n_freeR -= 1

        # create bond
        self.bp[mL][sL] = (mR, sR)
        self.bp[mR][sR] = (mL, sL)
        self.bond_list.append((mL, sL, mR, sR))
        self.total_bonds += 1

        # merge complexes if needed
        if not same:
            if len(self.cx_mols[cxL]) < len(self.cx_mols[cxR]):
                cxL, cxR = cxR, cxL
            for m in self.cx_mols[cxR]:
                self.cx_of[m] = cxL
            self.cx_freeL[cxL] += self.cx_freeL[cxR]
            self.cx_freeR[cxL] += self.cx_freeR[cxR]
            self.cx_mols[cxL] |= self.cx_mols[cxR]
            del self.cx_mols[cxR]
            del self.cx_freeL[cxR]
            del self.cx_freeR[cxR]

        cx = self.cx_of[mL]
        self.cx_freeL[cx] -= 1
        self.cx_freeR[cx] -= 1

        new_i = self.cx_freeL[cx] * self.cx_freeR[cx]
        self.intra_pairs += new_i - old_i

    def _remove_bond(self, idx):
        mL, sL, mR, sR = self.bond_list[idx]
        cx = self.cx_of[mL]

        # remove bond
        self.bp[mL][sL] = None
        self.bp[mR][sR] = None
        last = len(self.bond_list) - 1
        if idx != last:
            self.bond_list[idx] = self.bond_list[last]
        self.bond_list.pop()
        self.total_bonds -= 1

        # update molecule counts
        ob = self.nbonds[mL]
        self.nbonds[mL] -= 1
        self.nbonds[mR] -= 1

        free_after = 3 - self.nbonds[mL]
        if ob == 1:
            self.n_freeL_L1 -= (free_after - 1)
            self.n_freeL_L0 += free_after
        else:
            self.n_freeL_L1 += 1
        self.n_freeR += 1

        old_i = self.cx_freeL[cx] * self.cx_freeR[cx]
        self.cx_freeL[cx] += 1
        self.cx_freeR[cx] += 1

        # connectivity check: BFS from mL, does it reach mR?
        split = not self._reachable(mL, mR)

        if split:
            comp = self._component(mR)
            ncx = self.next_cx
            self.next_cx += 1
            fL = fR = 0
            for m in comp:
                self.cx_of[m] = ncx
                if self.mol_type[m] == 0:
                    fL += 3 - self.nbonds[m]
                else:
                    fR += 2 - self.nbonds[m]
            self.cx_mols[ncx] = comp
            self.cx_mols[cx] -= comp
            self.cx_freeL[cx] -= fL
            self.cx_freeR[cx] -= fR
            self.cx_freeL[ncx] = fL
            self.cx_freeR[ncx] = fR
            ni = (self.cx_freeL[cx] * self.cx_freeR[cx]
                  + self.cx_freeL[ncx] * self.cx_freeR[ncx])
        else:
            ni = self.cx_freeL[cx] * self.cx_freeR[cx]

        self.intra_pairs += ni - old_i

    def _reachable(self, src, dst):
        visited = set()
        visited.add(src)
        stack = [src]
        while stack:
            cur = stack.pop()
            for partner in self.bp[cur]:
                if partner is None:
                    continue
                nb = partner[0]
                if nb == dst:
                    return True
                if nb not in visited:
                    visited.add(nb)
                    stack.append(nb)
        return False

    def _component(self, start):
        visited = {start}
        stack = [start]
        while stack:
            cur = stack.pop()
            for partner in self.bp[cur]:
                if partner is None:
                    continue
                nb = partner[0]
                if nb not in visited:
                    visited.add(nb)
                    stack.append(nb)
        return visited

    # ── propensity & firing ───────────────────────────────────────────────

    def _propensities(self):
        # L0 molecules are unbound singletons → never in same cx as R
        p_cap   = KP1 * self.n_freeL_L0 * self.n_freeR
        inter1  = self.n_freeL_L1 * self.n_freeR - self.intra_pairs
        p_xlink = KP2 * max(inter1, 0)
        p_ring  = KP3 * self.intra_pairs
        p_unbind = KM * self.total_bonds
        return p_cap, p_xlink, p_ring, p_unbind

    def _fire(self, p_cap, p_xlink, p_ring, p_unbind):
        total = p_cap + p_xlink + p_ring + p_unbind
        r = random.random() * total

        if r < p_cap:
            self._fire_inter(l0_only=True)
        elif r < p_cap + p_xlink:
            self._fire_inter(l0_only=False)
        elif r < p_cap + p_xlink + p_ring:
            self._fire_ring()
        else:
            self._fire_unbind()

    def _pick_free_L_site(self, l0_only):
        """Pick a random free r-site on an L with the right bond count."""
        target_zero = l0_only
        candidates = []
        for m in range(LIG_TOT):
            if target_zero and self.nbonds[m] != 0:
                continue
            if not target_zero and self.nbonds[m] == 0:
                continue
            for s in range(3):
                if self.bp[m][s] is None:
                    candidates.append((m, s))
        return random.choice(candidates) if candidates else None

    def _pick_free_R_site(self):
        candidates = []
        for m in range(LIG_TOT, N_MOLS):
            for s in range(2):
                if self.bp[m][s] is None:
                    candidates.append((m, s))
        return random.choice(candidates) if candidates else None

    def _fire_inter(self, l0_only):
        # Rejection sampling: pick BOTH L and R fresh each attempt (unbiased)
        for _ in range(200):
            Ls = self._pick_free_L_site(l0_only)
            Rs = self._pick_free_R_site()
            if Ls and Rs and self.cx_of[Ls[0]] != self.cx_of[Rs[0]]:
                self._add_bond(Ls[0], Ls[1], Rs[0], Rs[1])
                return

    def _fire_ring(self):
        # weighted selection of a complex, then random pair within it
        candidates = []
        for cx_id in list(self.cx_mols):
            w = self.cx_freeL.get(cx_id, 0) * self.cx_freeR.get(cx_id, 0)
            if w > 0:
                candidates.append((cx_id, w))
        if not candidates:
            return
        total_w = sum(w for _, w in candidates)
        r = random.random() * total_w
        cum = 0
        sel_cx = candidates[0][0]
        for cx_id, w in candidates:
            cum += w
            if r < cum:
                sel_cx = cx_id
                break

        Ls, Rs = [], []
        for m in self.cx_mols[sel_cx]:
            if self.mol_type[m] == 0:
                for s in range(3):
                    if self.bp[m][s] is None:
                        Ls.append((m, s))
            else:
                for s in range(2):
                    if self.bp[m][s] is None:
                        Rs.append((m, s))
        if Ls and Rs:
            l = random.choice(Ls)
            r = random.choice(Rs)
            self._add_bond(l[0], l[1], r[0], r[1])

    def _fire_unbind(self):
        if self.bond_list:
            self._remove_bond(random.randrange(len(self.bond_list)))

    # ── observables ───────────────────────────────────────────────────────

    def observe(self):
        row = []

        # Ligand observables are Species type: count complexes containing
        # at least one L with the given bond count (not molecules).
        lig_cx = [set(), set(), set(), set()]   # bond_count -> set of cx_ids
        for m in range(LIG_TOT):
            lig_cx[self.nbonds[m]].add(self.cx_of[m])
        row.extend(len(s) for s in lig_cx)

        # Size observables: count complexes with exactly N R molecules
        size = [0] * (REC_TOT + 1)
        for cx_id, mols in self.cx_mols.items():
            nR = sum(1 for m in mols if self.mol_type[m] == 1)
            if 1 <= nR <= REC_TOT:
                size[nR] += 1
        row.extend(size[1:])
        return row

    # ── main loop ─────────────────────────────────────────────────────────

    def run(self):
        dt = T_END / N_STEPS
        results = [self.observe()]   # t=0

        next_out = dt
        out_idx = 1

        while out_idx <= N_STEPS:
            ps = self._propensities()
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

            self._fire(*ps)

        while len(results) <= N_STEPS:
            results.append(self.observe())

        return results


# ── Ensemble generation ──────────────────────────────────────────────────────

def run_ensemble(n_reps, base_seed=None):
    """Run n_reps simulations, return per-time-point arrays."""
    import numpy as np

    n_obs = len(OBS_NAMES)
    all_data = np.zeros((n_reps, N_STEPS + 1, n_obs))

    for rep in range(n_reps):
        seed = (base_seed + rep) if base_seed is not None else None
        t0 = timepkg.time()
        sim = Simulator(seed=seed)
        rows = sim.run()
        elapsed = timepkg.time() - t0
        print(f"  rep {rep+1}/{n_reps}  events≈{sim.total_bonds}  "
              f"bonds_final={sim.total_bonds}  {elapsed:.1f}s",
              file=sys.stderr)
        for ti, row in enumerate(rows):
            all_data[rep, ti, :] = row

    return all_data


def write_ensemble(all_data, out_dir):
    """Write mean.tsv, std.tsv, tint.tsv in the same format as NFsim ensemble."""
    import numpy as np
    os.makedirs(out_dir, exist_ok=True)

    n_reps = all_data.shape[0]
    dt = T_END / N_STEPS
    times = [dt * i for i in range(N_STEPS + 1)]

    mean = np.mean(all_data, axis=0)
    std  = np.std(all_data, axis=0, ddof=1)

    # mean.tsv
    with open(os.path.join(out_dir, "rm_tlbr_rings.mean.tsv"), "w") as f:
        f.write("time\t" + "\t".join(OBS_NAMES) + "\n")
        for ti in range(N_STEPS + 1):
            vals = "\t".join(f"{mean[ti, oi]:.6g}" for oi in range(len(OBS_NAMES)))
            f.write(f"{times[ti]:.6g}\t{vals}\n")

    # std.tsv
    with open(os.path.join(out_dir, "rm_tlbr_rings.std.tsv"), "w") as f:
        f.write("time\t" + "\t".join(OBS_NAMES) + "\n")
        for ti in range(N_STEPS + 1):
            vals = "\t".join(f"{std[ti, oi]:.6g}" for oi in range(len(OBS_NAMES)))
            f.write(f"{times[ti]:.6g}\t{vals}\n")

    # tint.tsv  (time-integrated mean and std, using trapezoidal rule)
    with open(os.path.join(out_dir, "rm_tlbr_rings.tint.tsv"), "w") as f:
        f.write("observable\tI_mean\tI_std\tn_reps\n")
        for oi, name in enumerate(OBS_NAMES):
            # per-rep time integrals (trapezoidal)
            integrals = np.zeros(n_reps)
            for r in range(n_reps):
                for ti in range(N_STEPS):
                    integrals[r] += 0.5 * (all_data[r, ti, oi] + all_data[r, ti+1, oi]) * dt
            f.write(f"{name}\t{np.mean(integrals):.6g}\t{np.std(integrals, ddof=1):.6g}"
                    f"\t{n_reps}\n")

    # Per-rep gdat files (needed by compute_noise_floor.py)
    rep_dir = os.path.join(
        os.path.dirname(out_dir), "replicates", "rm_tlbr_rings")
    os.makedirs(rep_dir, exist_ok=True)
    for r in range(n_reps):
        tag = f"rep_{r:03d}"
        gdat = os.path.join(rep_dir, f"{tag}_rm_tlbr_rings.gdat")
        with open(gdat, "w") as f:
            hdr = "#" + " ".join(f"{c:>15s}" for c in ["time"] + list(OBS_NAMES))
            f.write(hdr + "\n")
            for ti in range(N_STEPS + 1):
                vals = [f"{times[ti]:.8e}"] + [
                    f"{all_data[r, ti, oi]:.8e}" for oi in range(len(OBS_NAMES))]
                f.write("\t".join(vals) + "\n")
        # touch .done sentinel
        open(os.path.join(rep_dir, f"{tag}.done"), "w").close()

    print(f"Wrote ensemble to {out_dir}/", file=sys.stderr)
    print(f"Wrote {n_reps} per-rep gdat files to {rep_dir}/", file=sys.stderr)


def print_gdat(rows):
    """Print a single-rep gdat to stdout."""
    dt = T_END / N_STEPS
    header = "# time\t" + "\t".join(OBS_NAMES)
    print(header)
    for ti, row in enumerate(rows):
        t = dt * ti
        vals = "\t".join(str(v) for v in row)
        print(f"{t}\t{vals}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="rm_tlbr_rings Gillespie SSA")
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--output", type=str, default=None,
                    help="Output directory for ensemble files")
    ap.add_argument("--gdat", action="store_true",
                    help="Single rep, print gdat to stdout")
    args = ap.parse_args()

    if args.gdat:
        sim = Simulator(seed=args.seed)
        rows = sim.run()
        print_gdat(rows)
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
