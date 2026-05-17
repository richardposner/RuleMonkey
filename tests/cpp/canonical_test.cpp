// Direct unit tests for RuleMonkey's canonical complex labeler
// (cpp/rulemonkey/canonical.{hpp,cpp}; plan §6 steps 1-2).
//
// Two layers, per the plan's testing note (§6):
//
//   1. Hand-built ComplexGraph cases — pin specific shapes (single
//      molecule, asymmetric dimer, homo-dimer, within-molecule
//      symmetric components, chains, self-bond rings) to exact
//      canonical strings and to the input-order invariance the labeler
//      must hold.
//
//   2. A property-based test — the test that actually catches
//      canonicalization bugs.  It (a) generates a random connected
//      complex, applies a random graph isomorphism (permute molecules +
//      permute interchangeable same-name components), and asserts the
//      canonical label is unchanged; and (b) makes a targeted
//      structural change and asserts the label differs.
//
// Individualization-refinement (plan §3.2 step 3) is not yet
// implemented, so a genuinely symmetric complex reports
// `fully_refined == false`.  The property test asserts hard invariance
// only on fully-refined complexes and reports the residual count — when
// step 3 lands, every generated complex becomes fully refined and the
// reporting branch disappears.

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
// 1. Hand-built unit tests
// ===========================================================================

void test_single_molecule() {
  // A lone molecule with distinct, partly-stateful components renders
  // verbatim — components stay in declared slots, states with `~`.
  ComplexGraph g;
  g.add_molecule("A", comps({{"a", ""}, {"b", "x"}}));
  check_eq(canonical_label(g), "A(a,b~x)", "single molecule A(a,b~x)");
  check(canonicalize(g).fully_refined, "single molecule is fully refined");
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
  check(canonicalize(g1).fully_refined, "asymmetric dimer is fully refined");
}

void test_homodimer_symmetric() {
  // A(a!1).A(a!1): two molecules of one type, fully symmetric.  The
  // label is the same string regardless of build order even though the
  // complex is symmetric (so not fully_refined yet).
  ComplexGraph g;
  g.add_molecule("A", comps({{"a", ""}}));
  g.add_molecule("A", comps({{"a", ""}}));
  g.add_bond(0, 0, 1, 0);
  check_eq(canonical_label(g), "A(a!1).A(a!1)", "homodimer canonical string");
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
  check(canonicalize(g1).fully_refined,
        "A(a,a!1).B(b!1) is fully refined (WL distinguishes the two `a`s)");
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
  check(canonicalize(g1).fully_refined, "asymmetric 3-chain is fully refined");
}

void test_self_bond_ring() {
  // A single molecule whose two components are bonded to each other.
  ComplexGraph g;
  g.add_molecule("A", comps({{"l", ""}, {"r", ""}}));
  g.add_bond(0, 0, 0, 1);
  check_eq(canonical_label(g), "A(l!1,r!1)", "self-bonded molecule");
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
// 2. Property-based test
// ===========================================================================

// Molecule-type alphabet for the random generator.  Each entry is the
// ordered list of component names of that type; `A`/`B` have
// interchangeable components, `C` has distinct ones, `D` is monovalent.
// Held behind an accessor (function-local static) rather than a
// namespace-scope object so its allocating constructor never runs
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
// component.
bool generate(std::mt19937& rng, GenComplex& out) {
  out = GenComplex{};
  std::uniform_int_distribution<int> n_mol_d(1, 6);
  int const n_mol = n_mol_d(rng);

  std::uniform_int_distribution<int> type_d(0, static_cast<int>(types().size()) - 1);
  for (int m = 0; m < n_mol; ++m) {
    int const t = type_d(rng);
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

  // A few extra bonds add cycles (and self-bonds).
  std::uniform_int_distribution<int> extra_d(0, 2);
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
  int generated = 0, fully = 0, residual_match = 0, residual_total = 0;

  for (int run = 0; run < kRuns; ++run) {
    GenComplex gc;
    if (!generate(rng, gc))
      continue;
    ++generated;

    ComplexGraph g = build(gc);
    auto cf = canonicalize(g);

    // Determinism: the same input must always produce the same label.
    check_eq(canonical_label(g), cf.label, "canonicalize is deterministic");

    // (a) Isomorphism invariance.
    const GenComplex perm = permute(rng, gc);
    ComplexGraph gp = build(perm);
    auto cfp = canonicalize(gp);
    if (cf.fully_refined) {
      ++fully;
      check_eq(cfp.label, cf.label, "fully-refined canonical label is isomorphism-invariant");
      check(cfp.fully_refined, "fully_refined verdict is itself isomorphism-invariant");
    } else {
      // Symmetric residue — invariance is not yet guaranteed (step 3).
      ++residual_total;
      if (cfp.label == cf.label)
        ++residual_match;
    }

    // (b) A targeted structural change must change the label: flip the
    // state of one component.
    GenComplex tweaked = gc;
    // find any component and flip its state
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
               "property test: %d complexes generated, %d fully refined, "
               "%d symmetric residue (%d/%d permutation-stable anyway)\n",
               generated, fully, residual_total, residual_match, residual_total);
  check(generated > 1000, "generator produced a healthy sample");
  check(fully > 500, "a substantial fraction of complexes hit the fast path");
}

} // namespace

int main() {
  try {
    test_single_molecule();
    test_asymmetric_dimer_order_invariant();
    test_homodimer_symmetric();
    test_within_molecule_symmetric_components();
    test_chain_order_invariant();
    test_self_bond_ring();
    test_non_isomorphic_differ();
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
