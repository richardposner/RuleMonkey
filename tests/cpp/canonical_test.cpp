// Direct unit tests for RuleMonkey's canonical complex labeler
// (cpp/rulemonkey/canonical.{hpp,cpp}; plan §6 steps 1-3).
//
// Two layers, per the plan's testing note (§6):
//
//   1. Hand-built ComplexGraph cases — pin specific shapes (single
//      molecule, asymmetric dimer, chains, self-bond ring, and the
//      genuinely symmetric shapes that exercise individualization-
//      refinement: rings, homo-dimers/trimers, the ss_tlbr_rings and
//      ss_symmetric_homopoly feature-coverage topologies) to exact
//      canonical strings and to build-order invariance.
//
//   2. A property-based test — the test that actually catches
//      canonicalization bugs.  It (a) generates a random connected
//      complex, applies a random graph isomorphism (permute molecules +
//      permute interchangeable same-name components), and asserts the
//      canonical label is unchanged; and (b) makes a targeted
//      structural change and asserts the label differs.  The generator
//      is biased toward symmetric complexes (rings, homo-oligomers) so
//      the individualization-refinement search path is heavily
//      exercised.
//
// As of plan §3.2 step 3 the labeler is complete: individualization-
// refinement resolves genuine symmetry, so EVERY complex — symmetric or
// not — has a true canonical form.  The property test therefore asserts
// isomorphism invariance as a HARD invariant on every generated input
// (no fully-refined gate), and counts how many inputs needed the search.

#include "canonical.hpp"

#include <algorithm>
#include <array>
#include <cstdio>
#include <initializer_list>
#include <map>
#include <random>
#include <string>
#include <utility>
#include <vector>

using rulemonkey::canonical::canonical_label;
using rulemonkey::canonical::canonicalize;
using rulemonkey::canonical::ComplexGraph;

namespace {

int g_failures = 0;

void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

void check_eq(const std::string& got, const std::string& want, const std::string& msg) {
  if (got != want) {
    std::fprintf(stderr, "FAIL: %s\n  got:  %s\n  want: %s\n", msg.c_str(), got.c_str(),
                 want.c_str());
    ++g_failures;
  }
}

// Shorthand for an (name, state) component spec.
using C = std::pair<std::string, std::string>;
std::vector<C> comps(std::initializer_list<C> cs) { return {cs.begin(), cs.end()}; }

// ===========================================================================
// 1. Hand-built unit tests — asymmetric shapes (fast path)
// ===========================================================================

void test_single_molecule() {
  // A lone molecule with distinct, partly-stateful components renders
  // verbatim — components stay in declared slots, states with `~`.
  ComplexGraph g;
  g.add_molecule("A", comps({{"a", ""}, {"b", "x"}}));
  check_eq(canonical_label(g), "A(a,b~x)", "single molecule A(a,b~x)");
  check(canonicalize(g).fast_path, "single molecule takes the fast path");
}

void test_asymmetric_dimer_order_invariant() {
  // A(a!1).B(b!1): the two molecules have distinct types, so the
  // canonical order is fixed regardless of which is added first.
  ComplexGraph g1;
  g1.add_molecule("A", comps({{"a", ""}}));
  g1.add_molecule("B", comps({{"b", ""}}));
  g1.add_bond(0, 0, 1, 0);

  ComplexGraph g2; // built in the opposite molecule order
  g2.add_molecule("B", comps({{"b", ""}}));
  g2.add_molecule("A", comps({{"a", ""}}));
  g2.add_bond(0, 0, 1, 0);

  check_eq(canonical_label(g1), "A(a!1).B(b!1)", "asymmetric dimer canonical string");
  check_eq(canonical_label(g1), canonical_label(g2), "asymmetric dimer invariant to build order");
  check(canonicalize(g1).fast_path, "asymmetric dimer takes the fast path");
}

void test_within_molecule_symmetric_components() {
  // A has two interchangeable `a` components; one is bonded to B.
  // Which slot carries the bond must not change the canonical label.
  ComplexGraph g1;
  g1.add_molecule("A", comps({{"a", ""}, {"a", ""}}));
  g1.add_molecule("B", comps({{"b", ""}}));
  g1.add_bond(0, 0, 1, 0); // bond on A's slot 0

  ComplexGraph g2;
  g2.add_molecule("A", comps({{"a", ""}, {"a", ""}}));
  g2.add_molecule("B", comps({{"b", ""}}));
  g2.add_bond(0, 1, 1, 0); // bond on A's slot 1

  check_eq(canonical_label(g1), canonical_label(g2), "interchangeable component slot invariance");
  check(canonicalize(g1).fast_path, "A(a,a!1).B(b!1) takes the fast path (only one `a` is bonded)");
}

void test_chain_order_invariant() {
  // A 3-chain A(t!1).B(l!1,r!2).A(t!2): the two A's are distinguished
  // only through B, but B is asymmetric (distinct component names l/r),
  // so WL fully refines the whole complex.
  ComplexGraph g1;
  g1.add_molecule("A", comps({{"t", ""}}));
  g1.add_molecule("B", comps({{"l", ""}, {"r", ""}}));
  g1.add_molecule("A", comps({{"t", ""}}));
  g1.add_bond(0, 0, 1, 0);
  g1.add_bond(2, 0, 1, 1);

  ComplexGraph g2; // molecules added in a scrambled order
  g2.add_molecule("B", comps({{"l", ""}, {"r", ""}}));
  g2.add_molecule("A", comps({{"t", ""}}));
  g2.add_molecule("A", comps({{"t", ""}}));
  g2.add_bond(0, 0, 1, 0);
  g2.add_bond(0, 1, 2, 0);

  check_eq(canonical_label(g1), canonical_label(g2), "3-chain build-order invariance");
  check(canonicalize(g1).fast_path, "asymmetric 3-chain takes the fast path");
}

void test_self_bond_ring() {
  // A single molecule whose two distinctly-named components are bonded
  // to each other — asymmetric, fast path.
  ComplexGraph g;
  g.add_molecule("A", comps({{"l", ""}, {"r", ""}}));
  g.add_bond(0, 0, 0, 1);
  check_eq(canonical_label(g), "A(l!1,r!1)", "self-bonded molecule");
  check(canonicalize(g).fast_path, "self-bonded molecule takes the fast path");
}

void test_non_isomorphic_differ() {
  // The same shape, one component state flipped — must get a different
  // label so the species do not collapse in the dedup map.
  ComplexGraph g1;
  g1.add_molecule("A", comps({{"a", ""}}));
  g1.add_molecule("B", comps({{"b", ""}}));
  g1.add_bond(0, 0, 1, 0);

  ComplexGraph g2;
  g2.add_molecule("A", comps({{"a", "x"}}));
  g2.add_molecule("B", comps({{"b", ""}}));
  g2.add_bond(0, 0, 1, 0);

  check(canonical_label(g1) != canonical_label(g2), "state difference yields a distinct label");
}

// ===========================================================================
// 1b. Hand-built unit tests — symmetric shapes (individualization search)
// ===========================================================================

void test_homodimer_symmetric() {
  // A(a!1).A(a!1): two molecules of one type, fully symmetric.  1-WL
  // cannot break the molecule-vertex tie — individualization must.
  ComplexGraph g;
  g.add_molecule("A", comps({{"a", ""}}));
  g.add_molecule("A", comps({{"a", ""}}));
  g.add_bond(0, 0, 1, 0);
  check_eq(canonical_label(g), "A(a!1).A(a!1)", "homodimer canonical string");
  check(!canonicalize(g).fast_path, "homodimer needs the individualization search");
}

// Build a homo-ring of `n` copies of A(l,r): A[k].r bonds A[k+1].l, the
// last wraps to A[0].  `start` rotates the build so isomorphism
// invariance can be checked; `reverse` flips bond direction.
ComplexGraph make_lr_ring(int n, int start, bool reverse) {
  ComplexGraph g;
  for (int k = 0; k < n; ++k)
    g.add_molecule("A", comps({{"l", ""}, {"r", ""}}));
  for (int k = 0; k < n; ++k) {
    int const a = (start + k) % n;
    int const b = (start + k + 1) % n;
    if (reverse)
      g.add_bond(b, 0, a, 1); // b.l -- a.r
    else
      g.add_bond(a, 1, b, 0); // a.r -- b.l
  }
  return g;
}

void test_lr_ring_3() {
  // 3-membered ring of A(l,r) — a regular graph 1-WL cannot resolve.
  std::string const canon = canonical_label(make_lr_ring(3, 0, false));
  check_eq(canon, "A(l!1,r!2).A(l!2,r!3).A(l!3,r!1)", "3-ring canonical string");
  check(!canonicalize(make_lr_ring(3, 0, false)).fast_path, "3-ring needs the search");
  // Isomorphism invariance: every rotation and the reflection of the
  // ring must canonicalize to the same string.
  for (int s = 0; s < 3; ++s) {
    check_eq(canonical_label(make_lr_ring(3, s, false)), canon, "3-ring rotation invariance");
    check_eq(canonical_label(make_lr_ring(3, s, true)), canon, "3-ring reflection invariance");
  }
}

void test_lr_ring_4() {
  // 4-membered ring of A(l,r).  Molecules are emitted in refined-color
  // order (not ring-walk order), so the bond labels are not sequential.
  std::string const canon = canonical_label(make_lr_ring(4, 0, false));
  check_eq(canon, "A(l!1,r!2).A(l!3,r!1).A(l!2,r!4).A(l!4,r!3)", "4-ring canonical string");
  check(!canonicalize(make_lr_ring(4, 0, false)).fast_path, "4-ring needs the search");
  for (int s = 0; s < 4; ++s) {
    check_eq(canonical_label(make_lr_ring(4, s, false)), canon, "4-ring rotation invariance");
    check_eq(canonical_label(make_lr_ring(4, s, true)), canon, "4-ring reflection invariance");
  }
}

// ss_symmetric_homopoly shape: P(s,s) — one molecule type, two
// interchangeable `s` sites.  Build a chain of `n` copies, P[k].s bonds
// P[k+1].s; `ring` wraps the last back to P[0].  `start`/`reverse`
// rotate the build for invariance checks.
ComplexGraph make_homopoly(int n, bool ring, int start, bool reverse) {
  ComplexGraph g;
  for (int k = 0; k < n; ++k)
    g.add_molecule("P", comps({{"s", ""}, {"s", ""}}));
  int const links = ring ? n : n - 1;
  for (int e = 0; e < links; ++e) {
    int const a = (start + e) % n;
    int const b = (start + e + 1) % n;
    // s-slot 1 of `a` bonds s-slot 0 of `b` (slot choice is itself a
    // symmetry the canonicalizer must absorb).
    if (reverse)
      g.add_bond(b, 0, a, 1);
    else
      g.add_bond(a, 1, b, 0);
  }
  return g;
}

void test_homopoly_trimer_chain() {
  // ss_symmetric_homopoly: a linear homo-trimer P(s,s) — reflection-
  // symmetric (the two ends are interchangeable).  The doubly-bonded
  // center sorts first (no free site -> smaller molecule color).
  std::string const canon = canonical_label(make_homopoly(3, false, 0, false));
  check_eq(canon, "P(s!1,s!2).P(s,s!1).P(s,s!2)", "homopoly trimer chain canonical string");
  check(!canonicalize(make_homopoly(3, false, 0, false)).fast_path,
        "homopoly trimer chain needs the search");
  check_eq(canonical_label(make_homopoly(3, false, 0, true)), canon,
           "homopoly trimer chain reflection invariance");
}

void test_homopoly_ring_3() {
  // A 3-membered homo-ring of P(s,s) — every molecule and every site
  // interchangeable; the hardest small symmetric case.
  std::string const canon = canonical_label(make_homopoly(3, true, 0, false));
  check_eq(canon, "P(s!1,s!2).P(s!3,s!1).P(s!2,s!3)", "homopoly 3-ring canonical string");
  check(!canonicalize(make_homopoly(3, true, 0, false)).fast_path,
        "homopoly 3-ring needs the search");
  for (int s = 0; s < 3; ++s) {
    check_eq(canonical_label(make_homopoly(3, true, s, false)), canon,
             "homopoly 3-ring rotation invariance");
    check_eq(canonical_label(make_homopoly(3, true, s, true)), canon,
             "homopoly 3-ring reflection invariance");
  }
}

// ss_tlbr_rings shape: trivalent L(r,r,r) cross-linked with bivalent
// R(l,l).  The smallest ring alternates 2 L and 2 R; each L keeps one
// free `r`.  Modeled as a 4-cycle of ring positions 0..3 (even = L,
// odd = R); molecule index k is ring position (k+rot)%4, so `rot`
// rotates the build order — every rotation is the same ring.
ComplexGraph make_tlbr_ring(int rot) {
  ComplexGraph g;
  auto mol_of = [&](int p) { return (p - rot + 4) % 4; };
  for (int k = 0; k < 4; ++k) {
    int const p = (k + rot) % 4;
    if (p % 2 == 0)
      g.add_molecule("L", comps({{"r", ""}, {"r", ""}, {"r", ""}}));
    else
      g.add_molecule("R", comps({{"l", ""}, {"l", ""}}));
  }
  // Ring edge e joins ring positions e and e+1.  An L at position p
  // owns edges p-1 (its site 0) and p (its site 1); likewise for R.
  for (int e = 0; e < 4; ++e) {
    int const pa = e, pb = (e + 1) % 4;
    int const p_l = (pa % 2 == 0) ? pa : pb;
    int const p_r = (pa % 2 == 0) ? pb : pa;
    int const l_site = (e == p_l) ? 1 : 0;
    int const r_site = (e == p_r) ? 1 : 0;
    g.add_bond(mol_of(p_l), l_site, mol_of(p_r), r_site);
  }
  return g;
}

void test_tlbr_ring() {
  // 2L + 2R alternating ring — the ss_tlbr_rings minimal closure.
  std::string const canon = canonical_label(make_tlbr_ring(0));
  check_eq(canon, "L(r,r!1,r!2).R(l!1,l!3).R(l!4,l!2).L(r,r!3,r!4)", "tlbr ring canonical string");
  check(!canonicalize(make_tlbr_ring(0)).fast_path, "tlbr ring needs the search");
  // Rotating the build order is an isomorphism — same ring.
  for (int r = 0; r < 4; ++r)
    check_eq(canonical_label(make_tlbr_ring(r)), canon, "tlbr ring build-rotation invariance");
}

// ===========================================================================
// 2. Property-based test
// ===========================================================================

// Molecule-type alphabet for the random generator.  `A`/`B` have
// interchangeable components (the symmetry-prone types), `C` has
// distinct ones, `D` is monovalent.  Held behind an accessor
// (function-local static) so its allocating constructor never runs
// before main.
using TypeAlphabet = std::vector<std::pair<std::string, std::vector<std::string>>>;
const TypeAlphabet& types() {
  static const TypeAlphabet t = {
      {"A", {"s", "s", "s"}},
      {"B", {"s", "s"}},
      {"C", {"x", "y"}},
      {"D", {"s"}},
  };
  return t;
}

// Weighted type-index pool: heavily favors the interchangeable-component
// types A and B so generated complexes are often genuinely symmetric
// (rings, homo-oligomers) and the individualization search is the path
// being exercised, not the fast path.
const std::vector<int>& type_pool() {
  static const std::vector<int> pool = {0, 0, 0, 0, 1, 1, 1, 2, 3}; // A x4, B x3, C, D
  return pool;
}

// A generated complex, kept in a form the isomorphism applicator can
// permute: per-molecule type index + per-component state, plus a bond
// list addressed as (molecule, local component).
struct GenComplex {
  std::vector<int> type_idx;                    // per molecule
  std::vector<std::vector<std::string>> states; // [mol][local component]
  std::vector<std::array<int, 4>> bonds;        // {mol_a, local_a, mol_b, local_b}
};

ComplexGraph build(const GenComplex& gc) {
  ComplexGraph g;
  for (size_t m = 0; m < gc.type_idx.size(); ++m) {
    const auto& names = types()[gc.type_idx[m]].second;
    std::vector<C> cs;
    cs.reserve(names.size());
    for (size_t i = 0; i < names.size(); ++i)
      cs.emplace_back(names[i], gc.states[m][i]);
    g.add_molecule(types()[gc.type_idx[m]].first, cs);
  }
  for (const auto& b : gc.bonds)
    g.add_bond(b[0], b[1], b[2], b[3]);
  return g;
}

// Generate a random connected complex.  Returns false (caller retries)
// if a random spanning tree could not be completed for want of a free
// component.  Biased toward symmetric complexes: every so often all
// molecules share one type (homo-oligomers), and several extra bonds
// are added to close rings.
bool generate(std::mt19937& rng, GenComplex& out) {
  out = GenComplex{};
  std::uniform_int_distribution<int> n_mol_d(1, 7);
  int const n_mol = n_mol_d(rng);

  std::uniform_int_distribution<int> pool_d(0, static_cast<int>(type_pool().size()) - 1);
  // One run in three is a pure homo-oligomer: all molecules one type.
  std::uniform_int_distribution<int> mono_d(0, 2);
  bool const mono = mono_d(rng) == 0;
  int const mono_type = type_pool()[pool_d(rng)];

  for (int m = 0; m < n_mol; ++m) {
    int const t = mono ? mono_type : type_pool()[pool_d(rng)];
    out.type_idx.push_back(t);
    std::vector<std::string> st;
    std::uniform_int_distribution<int> state_d(0, 2); // 0 = stateless
    for (size_t i = 0; i < types()[t].second.size(); ++i) {
      int const s = state_d(rng);
      st.emplace_back(s == 0 ? "" : (s == 1 ? "p" : "q"));
    }
    out.states.push_back(std::move(st));
  }

  // free[m] = local component indices of molecule m not yet bonded.
  std::vector<std::vector<int>> free(n_mol);
  for (int m = 0; m < n_mol; ++m)
    for (size_t i = 0; i < types()[out.type_idx[m]].second.size(); ++i)
      free[m].push_back(static_cast<int>(i));

  auto take_free = [&](int m) -> int {
    std::uniform_int_distribution<int> d(0, static_cast<int>(free[m].size()) - 1);
    int const idx = d(rng);
    int const comp = free[m][idx];
    free[m].erase(free[m].begin() + idx);
    return comp;
  };

  // Spanning tree: connect molecule k to a random earlier molecule.
  for (int k = 1; k < n_mol; ++k) {
    std::vector<int> parents;
    for (int j = 0; j < k; ++j)
      if (!free[j].empty())
        parents.push_back(j);
    if (parents.empty() || free[k].empty())
      return false; // cannot keep the complex connected — retry
    std::uniform_int_distribution<int> pd(0, static_cast<int>(parents.size()) - 1);
    int const j = parents[pd(rng)];
    out.bonds.push_back({k, take_free(k), j, take_free(j)});
  }

  // Several extra bonds add cycles (and self-bonds) — biases the sample
  // toward ring-bearing, symmetric complexes.
  std::uniform_int_distribution<int> extra_d(0, 5);
  int const extra = extra_d(rng);
  for (int e = 0; e < extra; ++e) {
    std::vector<int> have;
    for (int m = 0; m < n_mol; ++m)
      if (!free[m].empty())
        have.push_back(m);
    if (have.empty())
      break;
    std::uniform_int_distribution<int> hd(0, static_cast<int>(have.size()) - 1);
    int const ma = have[hd(rng)];
    int const ca = take_free(ma);
    // pick a second endpoint (possibly the same molecule)
    have.clear();
    for (int m = 0; m < n_mol; ++m)
      if (!free[m].empty())
        have.push_back(m);
    if (have.empty())
      break;
    std::uniform_int_distribution<int> hd2(0, static_cast<int>(have.size()) - 1);
    int const mb = have[hd2(rng)];
    int const cb = take_free(mb);
    out.bonds.push_back({ma, ca, mb, cb});
  }
  return true;
}

// Apply a random graph isomorphism: permute molecules, and within each
// molecule permute components that share a name.  The result is the
// same species, so its canonical label must be unchanged.
GenComplex permute(std::mt19937& rng, const GenComplex& in) {
  int const n = static_cast<int>(in.type_idx.size());

  // pi: new molecule k <- old molecule pi[k].
  std::vector<int> pi(n);
  for (int i = 0; i < n; ++i)
    pi[i] = i;
  std::shuffle(pi.begin(), pi.end(), rng);
  std::vector<int> pi_inv(n);
  for (int k = 0; k < n; ++k)
    pi_inv[pi[k]] = k;

  // sigma[o]: name-preserving permutation of old molecule o's slots.
  // sigma[o][newlocal] = old local slot placed at newlocal.
  std::vector<std::vector<int>> sigma(n), sigma_inv(n);
  for (int o = 0; o < n; ++o) {
    const auto& names = types()[in.type_idx[o]].second;
    int const nc = static_cast<int>(names.size());
    std::vector<int> sig(nc);
    for (int i = 0; i < nc; ++i)
      sig[i] = i;
    // Shuffle each same-name group of slot positions among themselves.
    std::map<std::string, std::vector<int>> groups;
    for (int i = 0; i < nc; ++i)
      groups[names[i]].push_back(i);
    for (auto& [name, slots] : groups) {
      std::vector<int> shuffled = slots;
      std::shuffle(shuffled.begin(), shuffled.end(), rng);
      for (size_t i = 0; i < slots.size(); ++i)
        sig[slots[i]] = shuffled[i]; // newlocal slots[i] <- old slot shuffled[i]
    }
    sigma[o] = sig;
    std::vector<int> inv(nc);
    for (int i = 0; i < nc; ++i)
      inv[sig[i]] = i;
    sigma_inv[o] = inv;
  }

  GenComplex out;
  out.type_idx.resize(n);
  out.states.resize(n);
  for (int k = 0; k < n; ++k) {
    int const o = pi[k];
    out.type_idx[k] = in.type_idx[o];
    int const nc = static_cast<int>(in.states[o].size());
    out.states[k].resize(nc);
    for (int l = 0; l < nc; ++l)
      out.states[k][l] = in.states[o][sigma[o][l]];
  }
  for (const auto& b : in.bonds) {
    // old (mol o, local p) -> new (mol pi_inv[o], local sigma_inv[o][p])
    int const oa = b[0], ob = b[2];
    out.bonds.push_back({pi_inv[oa], sigma_inv[oa][b[1]], pi_inv[ob], sigma_inv[ob][b[3]]});
  }
  return out;
}

void test_property_based() {
  // Fixed seed: a property test must be reproducible so a failure is
  // re-runnable. NOLINTNEXTLINE(bugprone-random-generator-seed)
  std::mt19937 rng(0xC0FFEEU);
  int const kRuns = 4000;
  int generated = 0, fast = 0, searched = 0;

  for (int run = 0; run < kRuns; ++run) {
    GenComplex gc;
    if (!generate(rng, gc))
      continue;
    ++generated;

    ComplexGraph g = build(gc);
    auto cf = canonicalize(g);
    if (cf.fast_path)
      ++fast;
    else
      ++searched;

    // Determinism: the same input must always produce the same label.
    check_eq(canonical_label(g), cf.label, "canonicalize is deterministic");

    // (a) Isomorphism invariance — a HARD invariant on every input now
    // that individualization-refinement resolves all symmetry (step 3).
    const GenComplex perm = permute(rng, gc);
    ComplexGraph gp = build(perm);
    auto cfp = canonicalize(gp);
    check_eq(cfp.label, cf.label, "canonical label is isomorphism-invariant");
    check(cfp.fast_path == cf.fast_path, "fast_path verdict is itself isomorphism-invariant");

    // (b) A targeted structural change must change the label: flip the
    // state of one component.
    GenComplex tweaked = gc;
    {
      int const m = run % static_cast<int>(tweaked.states.size());
      int const c = 0;
      std::string& s = tweaked.states[m][c];
      s = (s.empty() ? "p" : (s == "p" ? "q" : ""));
    }
    ComplexGraph gt = build(tweaked);
    check(canonical_label(gt) != cf.label, "a flipped component state yields a distinct label");
  }

  std::fprintf(stderr,
               "property test: %d complexes generated, %d fast path, "
               "%d needed individualization search\n",
               generated, fast, searched);
  check(generated > 1000, "generator produced a healthy sample");
  check(fast > 200, "a substantial fraction of complexes hit the fast path");
  check(searched > 200, "the symmetric bias exercises the individualization search");
}

} // namespace

int main() {
  try {
    test_single_molecule();
    test_asymmetric_dimer_order_invariant();
    test_within_molecule_symmetric_components();
    test_chain_order_invariant();
    test_self_bond_ring();
    test_non_isomorphic_differ();
    test_homodimer_symmetric();
    test_lr_ring_3();
    test_lr_ring_4();
    test_homopoly_trimer_chain();
    test_homopoly_ring_3();
    test_tlbr_ring();
    test_property_based();
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: canonical labeler (hand-built shapes + isomorphism "
                       "property test) all pass\n");
  return 0;
}
