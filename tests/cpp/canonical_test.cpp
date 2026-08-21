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
//
// The integer entry point (GH #53) rides on the same two layers.  Every
// hand-built shape and every generated complex is ALSO put through
// `canonical_order_fast` over the RankedComplex the `rank_*` tables
// build from its own names, and the answer is pinned against
// `canonicalize().mol_order`.  That is the check that matters for #53:
// the two paths derive their initial colors completely differently —
// one interns strings, the other ranks integers — and the whole fix
// rests on those two rankings being the same ranking.

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
using rulemonkey::canonical::canonical_order;
using rulemonkey::canonical::canonical_order_fast;
using rulemonkey::canonical::canonicalize;
using rulemonkey::canonical::ComplexGraph;
using rulemonkey::canonical::rank_component_names;
using rulemonkey::canonical::rank_molecule_type_names;
using rulemonkey::canonical::rank_state_names;
using rulemonkey::canonical::RankedComplex;
using rulemonkey::canonical::Workspace;

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

// The `.`-separated molecule pieces of a canonical label.  Splitting on
// `.` is safe: a BNGL molecule never contains one.
std::vector<std::string> label_pieces(const std::string& label) {
  std::vector<std::string> out;
  size_t start = 0;
  while (true) {
    size_t const dot = label.find('.', start);
    out.push_back(label.substr(start, dot - start));
    if (dot == std::string::npos)
      break;
    start = dot + 1;
  }
  return out;
}

// ---------------------------------------------------------------------------
// The integer entry point, driven exactly as the engine drives it: rank
// the graph's own names through the canonicalizer's tables, restate the
// complex as a RankedComplex, ask canonical_order_fast.
//
// The rank_* tables take names in the caller's indexing and tolerate
// duplicates, so the graph's flat component array can be ranked whole.
// ---------------------------------------------------------------------------
RankedComplex to_ranked(const ComplexGraph& g) {
  std::vector<std::string> types, names, states;
  types.reserve(g.molecules.size());
  names.reserve(g.components.size());
  states.reserve(g.components.size());
  for (const auto& m : g.molecules)
    types.push_back(m.type_name);
  for (const auto& c : g.components) {
    names.push_back(c.name);
    states.push_back(c.state);
  }
  const auto type_rank = rank_molecule_type_names(types);
  const auto name_rank = rank_component_names(names);
  const auto state_rank = rank_state_names(states);

  RankedComplex rc;
  rc.molecules.resize(g.molecules.size());
  for (size_t m = 0; m < g.molecules.size(); ++m)
    rc.molecules[m] = {type_rank[m], g.molecules[m].first_comp, g.molecules[m].n_comp};
  rc.components.resize(g.components.size());
  for (size_t c = 0; c < g.components.size(); ++c)
    rc.components[c] = {name_rank[c], state_rank[c], g.components[c].partner};
  return rc;
}

// Whenever the integer path answers, it must answer exactly what
// canonicalize() does; and it must answer at least everywhere the
// string fast path does, since its gate is the weaker of the two.
// Returns whether it answered, so a caller can pin the fallback too.
bool check_ranked_agrees(const ComplexGraph& g, const std::string& msg) {
  static Workspace ws; // one workspace, reused — the shape #53 is about
  const auto cf = canonicalize(g);
  const RankedComplex rc = to_ranked(g);
  std::vector<int> order;
  bool const answered = canonical_order_fast(rc, ws, order);
  if (answered)
    check(order == cf.mol_order, "canonical_order_fast matches mol_order: " + msg);
  else
    check(!cf.fast_path, "canonical_order_fast declines only a searched complex: " + msg);
  return answered;
}

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

  // mol_order names the INPUT molecule written at each position of the
  // label, so it inverts with the build order while the label does not.
  check(canonicalize(g1).mol_order == std::vector<int>{0, 1}, "dimer mol_order, A added first");
  check(canonicalize(g2).mol_order == std::vector<int>{1, 0}, "dimer mol_order, B added first");
}

void test_mol_order_picks_the_same_subunit() {
  // GH #52's shape: a homodimer whose two subunits differ only in the
  // state of a free component.  Whichever order the pool happens to hold
  // them in, mol_order must name the SAME subunit first — that is what
  // lets a caller price one instance per complex without the answer
  // depending on how the complex was assembled.  BNG2 prices this
  // species at its `m~0` subunit; RuleMonkey's canonical order agrees,
  // and this pins that agreement.
  const auto build_dimer = [](bool mod_first) {
    ComplexGraph g;
    g.add_molecule("Wz", comps({{"d", ""}, {"m", mod_first ? "1" : "0"}}));
    g.add_molecule("Wz", comps({{"d", ""}, {"m", mod_first ? "0" : "1"}}));
    g.add_bond(0, 0, 1, 0);
    return g;
  };

  ComplexGraph const g_unmod_first = build_dimer(false);
  ComplexGraph const g_mod_first = build_dimer(true);
  check_eq(canonical_label(g_unmod_first), "Wz(d!1,m~0).Wz(d!1,m~1)", "mixed dimer canonical form");
  check_eq(canonical_label(g_mod_first), canonical_label(g_unmod_first),
           "mixed dimer label invariant to build order");

  // Position 0 is the `m~0` subunit in both, which is input molecule 0
  // in one build order and input molecule 1 in the other.
  check(canonicalize(g_unmod_first).mol_order.front() == 0, "m~0 subunit leads, unmod built first");
  check(canonicalize(g_mod_first).mol_order.front() == 1, "m~0 subunit leads, mod built first");
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

// A chain of `n` copies of A(l,r): position p's `r` bonds position p+1's
// `l`, and the chain does NOT close, so the two ends carry a free `l` and
// a free `r`.  `rot` shifts the position -> molecule-index map, so every
// rot builds the same chain in a different order.
ComplexGraph make_lr_chain(int n, int rot) {
  ComplexGraph g;
  auto mol_of = [&](int p) { return (p + rot) % n; };
  for (int k = 0; k < n; ++k)
    g.add_molecule("A", comps({{"l", ""}, {"r", ""}}));
  for (int p = 0; p + 1 < n; ++p)
    g.add_bond(mol_of(p), 1, mol_of(p + 1), 0); // p.r -- (p+1).l
  return g;
}

void test_lr_chain_refines_all_the_way() {
  // The shape that asks the most of refinement.  `l` and `r` are different
  // component names, so reversing the chain is not an automorphism and the
  // complex is asymmetric — but 1-WL only learns that one hop per round, so
  // the colors have to walk in from the two ends and a chain of n needs on
  // the order of n rounds to separate every molecule.
  //
  // That is what makes it the regression test for GH #56's worklist
  // refinement, which skips the cells a round cannot split: a cell it
  // wrongly fails to mark simply stops splitting, and refinement then
  // returns a partition that is too coarse.  Nothing about that is loud —
  // the individualization search still finds a correct canonical form — so
  // the assertion that catches it is `fast_path`, which says 1-WL alone
  // separated the complex and is FALSE the moment refinement gives up
  // early.  A 24-chain needs ~24 rounds, deep enough that a worklist that
  // ran dry after a couple would be caught here.
  // Molecules are emitted in refined-color order, not chain-walk order, so
  // the bond labels come out shuffled — same as the 4-ring below.
  check_eq(canonical_label(make_lr_chain(4, 0)), "A(l!1,r!2).A(l!3,r!1).A(l,r!3).A(l!2,r)",
           "4-chain canonical string");
  for (int const n : {2, 3, 4, 8, 24}) {
    std::string const canon = canonical_label(make_lr_chain(n, 0));
    check(canonicalize(make_lr_chain(n, 0)).fast_path,
          "the " + std::to_string(n) + "-chain refines without the search");
    for (int rot = 0; rot < n; ++rot) {
      check_eq(canonical_label(make_lr_chain(n, rot)), canon,
               "the " + std::to_string(n) + "-chain is invariant to build order");
      check(check_ranked_agrees(make_lr_chain(n, rot),
                                std::to_string(n) + "-chain rot " + std::to_string(rot)),
            "the " + std::to_string(n) + "-chain answers on the int path");
    }
  }

  // The same chain carrying a per-molecule state, which is the live shape
  // (GH #52's tagged catalyst, edited by a toggle rule).  States break the
  // symmetry sooner in some places and later in others, so the refinement
  // ends up with cells splitting in different rounds — the case where a
  // stale worklist entry would matter.
  for (int pattern = 0; pattern < 8; ++pattern) {
    ComplexGraph g;
    int const n = 12;
    for (int k = 0; k < n; ++k)
      g.add_molecule("W",
                     comps({{"l", ""}, {"r", ""}, {"m", ((pattern >> (k % 3)) & 1) ? "1" : "0"}}));
    for (int p = 0; p + 1 < n; ++p)
      g.add_bond(p, 1, p + 1, 0);
    check(canonicalize(g).fast_path, "a state-carrying 12-chain refines without the search");
    check(check_ranked_agrees(g, "state-carrying 12-chain"),
          "a state-carrying 12-chain answers on the int path");
  }
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
// 1a. The integer entry point — the rank rules that could silently
//     diverge from the color strings (GH #53)
// ===========================================================================

void test_ranked_rank_rules() {
  // A component name is written `name~state`, and `~` outranks every
  // character a BNGL identifier may hold, so `ab~` sorts BEFORE `a~`.
  // Rank component names in plain lexicographic order instead and this
  // complex comes out in the opposite order — the two molecules differ
  // only in which of their prefix-named components carries the state.
  //
  //   M0 color: M:E(ab~,a~p)     M1 color: M:E(ab~p,a~)
  //             -> M0 < M1, because `,` < `p`
  //
  // Under a plain-lexicographic component-name rank the keys would read
  // (a,p),(ab,"") and (a,""),(ab,p), putting M1 first.  So this shape
  // fails loudly if rank_component_names ever loses its `~`.
  ComplexGraph g;
  g.add_molecule("E", comps({{"a", "p"}, {"ab", ""}, {"l", ""}}));
  g.add_molecule("E", comps({{"a", ""}, {"ab", "p"}, {"l", ""}}));
  g.add_bond(0, 2, 1, 2);
  check(canonicalize(g).mol_order == std::vector<int>{0, 1},
        "prefix component names order M0 first");
  check(check_ranked_agrees(g, "prefix component names"),
        "prefix-name dimer answers on the int path");

  // The same shape with the state carried by the *shorter*-named
  // component on the other molecule — the mirror image, to catch a rank
  // rule that happens to agree on one arrangement by luck.
  ComplexGraph gm;
  gm.add_molecule("E", comps({{"a", ""}, {"ab", "p"}, {"l", ""}}));
  gm.add_molecule("E", comps({{"a", "p"}, {"ab", ""}, {"l", ""}}));
  gm.add_bond(0, 2, 1, 2);
  check(canonicalize(gm).mol_order == std::vector<int>{1, 0}, "prefix component names, mirrored");
  check_ranked_agrees(gm, "prefix component names, mirrored");

  // A molecule type name is written `Type(`, and `(` is outranked by
  // every identifier character, so a prefix type name sorts FIRST —
  // the opposite convention from component names.
  ComplexGraph gt;
  gt.add_molecule("Ez", comps({{"l", ""}}));
  gt.add_molecule("E", comps({{"a", ""}, {"ab", ""}, {"l", ""}}));
  gt.add_bond(0, 0, 1, 2);
  check(canonicalize(gt).mol_order == std::vector<int>{1, 0},
        "prefix molecule type names order E first");
  check_ranked_agrees(gt, "prefix molecule type names");

  // A state is written last, so a prefix state sorts first, and a
  // component with no internal state ("") sorts ahead of every state.
  ComplexGraph gs;
  gs.add_molecule("C", comps({{"x", "pq"}, {"y", ""}}));
  gs.add_molecule("C", comps({{"x", "p"}, {"y", ""}}));
  gs.add_bond(0, 1, 1, 1);
  check(canonicalize(gs).mol_order == std::vector<int>{1, 0}, "prefix states order `p` first");
  check_ranked_agrees(gs, "prefix states");

  // A shorter free-component list is a prefix of a longer one and sorts
  // first, exactly as `)` < `,` makes it in the color string.  F0 has
  // `b` taken by the D, so its color is `M:F(a~)` against F1's
  // `M:F(a~,b~)` and F0 must be written first.
  ComplexGraph gf;
  gf.add_molecule("F", comps({{"a", ""}, {"b", ""}, {"l", ""}}));
  gf.add_molecule("F", comps({{"a", ""}, {"b", ""}, {"l", ""}}));
  gf.add_molecule("D", comps({{"s", ""}}));
  gf.add_bond(0, 2, 1, 2); // the two F's are joined through `l`
  gf.add_bond(0, 1, 2, 0); // ...and F0's `b` is taken, shortening its free list
  const auto gf_order = canonicalize(gf).mol_order;
  const auto at = [&](int mol) {
    return std::find(gf_order.begin(), gf_order.end(), mol) - gf_order.begin();
  };
  check(at(0) < at(1), "the shorter free-component list leads");
  check_ranked_agrees(gf, "free-component list length");

  // Degenerate inputs the engine can hand it: an empty complex, and a
  // lone molecule (which the engine short-circuits, but the entry point
  // must still answer).
  ComplexGraph const empty;
  std::vector<int> order{7, 7, 7};
  Workspace ws;
  check(canonical_order_fast(to_ranked(empty), ws, order), "empty complex answers");
  check(order.empty(), "empty complex has an empty order");
  ComplexGraph one;
  one.add_molecule("A", comps({{"a", ""}, {"b", "x"}}));
  check(check_ranked_agrees(one, "single molecule"), "single molecule answers on the int path");
}

void test_ranked_declines_symmetric() {
  // A homodimer's two molecules are interchangeable, so refinement
  // leaves them sharing a color and only the RENDER can pick between
  // them.  The integer path has no strings and must say so.
  ComplexGraph g;
  g.add_molecule("A", comps({{"s", ""}, {"s", ""}, {"s", ""}}));
  g.add_molecule("A", comps({{"s", ""}, {"s", ""}, {"s", ""}}));
  g.add_bond(0, 0, 1, 0);
  Workspace ws;
  std::vector<int> order;
  check(!canonical_order_fast(to_ranked(g), ws, order), "homodimer declines on the int path");

  // But a molecule with two interchangeable bonded components is NOT a
  // reason to decline: the molecules already sit in distinct color
  // classes, so no leaf below can reorder them even though the string
  // path still has to search to pick a render.
  ComplexGraph gt;
  gt.add_molecule("A", comps({{"s", ""}, {"s", ""}, {"s", ""}}));
  gt.add_molecule("B", comps({{"s", ""}, {"s", ""}}));
  gt.add_bond(0, 0, 1, 0);
  gt.add_bond(0, 1, 1, 1);
  check(!canonicalize(gt).fast_path, "the two-bond A/B pair still searches for its label");
  check(check_ranked_agrees(gt, "two interchangeable bonds"),
        "a within-molecule tie alone does not make the int path decline");
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
// distinct ones, `D` is monovalent.  `E`/`Ez` carry the name shapes
// that separate the color-string order from a naive one (GH #53): a
// component name that is a prefix of another, and a type name that is a
// prefix of another.  Held behind an accessor (function-local static)
// so its allocating constructor never runs before main.
using TypeAlphabet = std::vector<std::pair<std::string, std::vector<std::string>>>;
const TypeAlphabet& types() {
  static const TypeAlphabet t = {
      {"A", {"s", "s", "s"}}, {"B", {"s", "s"}},  {"C", {"x", "y"}},
      {"D", {"s"}},           {"E", {"a", "ab"}}, {"Ez", {"a", "ab"}},
  };
  return t;
}

// Weighted type-index pool: heavily favors the interchangeable-component
// types A and B so generated complexes are often genuinely symmetric
// (rings, homo-oligomers) and the individualization search is the path
// being exercised, not the fast path.
const std::vector<int>& type_pool() {
  // A x6, B x5, C, D, E, Ez — A and B keep the same share of the pool
  // they had before E/Ez joined it, so the symmetric bias is unchanged.
  static const std::vector<int> pool = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 3, 4, 5};
  return pool;
}

// Uniform draw from [0, n), and a Fisher-Yates shuffle, both hand-rolled.
//
// std::uniform_int_distribution and std::shuffle are deterministic for a
// given seed but NOT portable: libc++ and libstdc++ turn the same engine
// output into different draws, so the same seed walks a different sample
// on macOS than on Linux.  That is fine for "a random complex" and fatal
// for anything that pins WHICH complexes the run produced — the fast-path
// / search split below is exactly such a pin, and it is what says
// refinement is as powerful as it was (GH #56).  These two make the whole
// property test reproducible across platforms instead of only across runs
// on one.  Modulo bias is irrelevant at these ranges.
int pick(std::mt19937& rng, int n) { return static_cast<int>(rng() % static_cast<unsigned>(n)); }

template <class T> void shuffle_v(std::vector<T>& v, std::mt19937& rng) {
  for (size_t i = v.size(); i > 1; --i)
    std::swap(v[i - 1], v[static_cast<size_t>(pick(rng, static_cast<int>(i)))]);
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
  int const n_mol = 1 + pick(rng, 7);

  int const n_pool = static_cast<int>(type_pool().size());
  // One run in three is a pure homo-oligomer: all molecules one type.
  bool const mono = pick(rng, 3) == 0;
  int const mono_type = type_pool()[pick(rng, n_pool)];

  for (int m = 0; m < n_mol; ++m) {
    int const t = mono ? mono_type : type_pool()[pick(rng, n_pool)];
    out.type_idx.push_back(t);
    std::vector<std::string> st;
    // 0 = stateless; `p`/`pq` are a prefix pair, so state ordering is
    // exercised at the one place it could go wrong (GH #53).
    for (size_t i = 0; i < types()[t].second.size(); ++i) {
      int const s = pick(rng, 4);
      st.emplace_back(s == 0 ? "" : (s == 1 ? "p" : (s == 2 ? "q" : "pq")));
    }
    out.states.push_back(std::move(st));
  }

  // free[m] = local component indices of molecule m not yet bonded.
  std::vector<std::vector<int>> free(n_mol);
  for (int m = 0; m < n_mol; ++m)
    for (size_t i = 0; i < types()[out.type_idx[m]].second.size(); ++i)
      free[m].push_back(static_cast<int>(i));

  auto take_free = [&](int m) -> int {
    int const idx = pick(rng, static_cast<int>(free[m].size()));
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
    int const j = parents[pick(rng, static_cast<int>(parents.size()))];
    out.bonds.push_back({k, take_free(k), j, take_free(j)});
  }

  // Several extra bonds add cycles (and self-bonds) — biases the sample
  // toward ring-bearing, symmetric complexes.
  int const extra = pick(rng, 6);
  for (int e = 0; e < extra; ++e) {
    std::vector<int> have;
    for (int m = 0; m < n_mol; ++m)
      if (!free[m].empty())
        have.push_back(m);
    if (have.empty())
      break;
    int const ma = have[pick(rng, static_cast<int>(have.size()))];
    int const ca = take_free(ma);
    // pick a second endpoint (possibly the same molecule)
    have.clear();
    for (int m = 0; m < n_mol; ++m)
      if (!free[m].empty())
        have.push_back(m);
    if (have.empty())
      break;
    int const mb = have[pick(rng, static_cast<int>(have.size()))];
    int const cb = take_free(mb);
    out.bonds.push_back({ma, ca, mb, cb});
  }
  return true;
}

// Apply a random graph isomorphism: permute molecules, and within each
// molecule permute components that share a name.  The result is the
// same species, so its canonical label must be unchanged.
// `pi_out`, when given, receives the molecule permutation: new molecule
// k is old molecule (*pi_out)[k].
GenComplex permute(std::mt19937& rng, const GenComplex& in, std::vector<int>* pi_out = nullptr) {
  int const n = static_cast<int>(in.type_idx.size());

  // pi: new molecule k <- old molecule pi[k].
  std::vector<int> pi(n);
  for (int i = 0; i < n; ++i)
    pi[i] = i;
  shuffle_v(pi, rng);
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
      shuffle_v(shuffled, rng);
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
  if (pi_out != nullptr)
    *pi_out = pi;
  return out;
}

void test_property_based() {
  // Fixed seed, and `pick`/`shuffle_v` rather than the standard
  // distributions, so the sample is the same everywhere and a failure
  // reported by one platform's CI is re-runnable on another's.
  // NOLINTNEXTLINE(bugprone-random-generator-seed)
  std::mt19937 rng(0xC0FFEEU);
  int const kRuns = 4000;
  int generated = 0, fast = 0, searched = 0;
  int ranked_answered = 0, ranked_declined = 0;

  for (int run = 0; run < kRuns; ++run) {
    GenComplex gc;
    if (!generate(rng, gc))
      continue;
    ++generated;

    ComplexGraph const g = build(gc);
    auto cf = canonicalize(g);
    if (cf.fast_path)
      ++fast;
    else
      ++searched;

    // Determinism: the same input must always produce the same label.
    check_eq(canonical_label(g), cf.label, "canonicalize is deterministic");

    // The order-only entry point skips the fast path's render.  It must
    // still answer identically — this is the check that keeps the two
    // from drifting apart, since only one of them builds the string.
    check(canonical_order(g) == cf.mol_order, "canonical_order matches canonicalize().mol_order");

    // Same question of the integer entry point, which reaches its
    // initial colors through the rank_* tables rather than through the
    // color strings (GH #53).  It may decline — a complex with
    // interchangeable molecules needs the render to choose — but when
    // it answers, it answers the same thing.
    if (check_ranked_agrees(g, "generated complex"))
      ++ranked_answered;
    else
      ++ranked_declined;

    // (a) Isomorphism invariance — a HARD invariant on every input now
    // that individualization-refinement resolves all symmetry (step 3).
    std::vector<int> pi;
    const GenComplex perm = permute(rng, gc, &pi);
    ComplexGraph const gp = build(perm);
    auto cfp = canonicalize(gp);
    check_eq(cfp.label, cf.label, "canonical label is isomorphism-invariant");
    check(cfp.fast_path == cf.fast_path, "fast_path verdict is itself isomorphism-invariant");
    if (check_ranked_agrees(gp, "permuted complex"))
      ++ranked_answered;
    else
      ++ranked_declined;

    // (a2) mol_order names the molecules the label writes, in order.
    // The label's k-th `.`-piece must be molecule mol_order[k]'s, so at
    // minimum its type name has to match — the check that catches an
    // inverted or off-by-one permutation.
    {
      const auto pieces = label_pieces(cf.label);
      bool const sized =
          cf.mol_order.size() == gc.type_idx.size() && pieces.size() == cf.mol_order.size();
      check(sized, "mol_order has one entry per molecule");
      if (sized) {
        std::vector<char> seen(cf.mol_order.size(), 0);
        for (size_t k = 0; k < cf.mol_order.size(); ++k) {
          int const m = cf.mol_order[k];
          bool const in_range = m >= 0 && m < static_cast<int>(seen.size());
          check(in_range && seen[m] == 0, "mol_order is a permutation of the molecules");
          if (!in_range)
            break;
          seen[m] = 1;
          const std::string& want = types()[gc.type_idx[m]].first;
          check(pieces[k].compare(0, want.size(), want) == 0 && pieces[k][want.size()] == '(',
                "mol_order[k] is the molecule the label writes at position k");
        }
      }
    }

    // (a3) On the fast path every molecule vertex ends in its own
    // refined color, so the canonical order is FORCED — there is no
    // automorphism left to permute it.  The two orderings must then
    // name corresponding molecules, element for element.  (Off the fast
    // path the complex has genuine symmetry and interchangeable
    // molecules may legitimately swap, so only (a2) applies.)
    if (cf.fast_path && cfp.mol_order.size() == cf.mol_order.size()) {
      bool matched = true;
      for (size_t k = 0; k < cf.mol_order.size(); ++k)
        matched = matched && pi[cfp.mol_order[k]] == cf.mol_order[k];
      check(matched, "fast-path mol_order tracks the isomorphism");
    }

    // (b) A targeted structural change must change the label: flip the
    // state of one component.
    GenComplex tweaked = gc;
    {
      int const m = run % static_cast<int>(tweaked.states.size());
      int const c = 0;
      std::string& s = tweaked.states[m][c];
      s = (s.empty() ? "p" : (s == "p" ? "q" : ""));
    }
    ComplexGraph const gt = build(tweaked);
    check(canonical_label(gt) != cf.label, "a flipped component state yields a distinct label");
  }

  std::fprintf(stderr,
               "property test: %d complexes generated, %d fast path, "
               "%d needed individualization search; integer entry point "
               "answered %d, declined %d\n",
               generated, fast, searched, ranked_answered, ranked_declined);
  check(generated > 1000, "generator produced a healthy sample");
  check(fast > 200, "a substantial fraction of complexes hit the fast path");
  check(searched > 200, "the symmetric bias exercises the individualization search");
  check(ranked_answered > 200, "the integer entry point answers on a healthy sample");
  // Declines are rare on purpose: the gate is "molecules in distinct
  // color classes", so a complex that needs the search to pick a LABEL
  // usually still has its ORDER fixed already.  The deterministic
  // homodimer in test_ranked_declines_symmetric is what pins the
  // fallback's behaviour; this only asserts the sample reaches it.
  check(ranked_declined > 20, "…and its decline path is exercised too");

  // The bounds above say the sample is healthy.  These say refinement is
  // exactly as POWERFUL as it was, which the bounds cannot: the generator
  // is seeded and canonicalization is deterministic, so the split between
  // the fast path and the search is a property of the code under test, and
  // refinement that stops one round early moves complexes across it
  // without failing anything else — every label stays a correct canonical
  // form, just a different one (GH #56).  The same goes for the integer
  // entry point, which declines exactly when refinement leaves two
  // molecules sharing a color.
  //
  // Re-record all four from the line the test prints if the generator or
  // its alphabet changes; do NOT re-record them for a change to
  // canonicalization itself without knowing which complexes moved and why.
  check(fast == 3586, "fast-path count is unchanged");
  check(searched == 345, "individualization-search count is unchanged");
  check(ranked_answered == 7788, "integer entry point answers as often as before");
  check(ranked_declined == 74, "…and declines as often as before");
}

} // namespace

int main() {
  try {
    test_single_molecule();
    test_asymmetric_dimer_order_invariant();
    test_mol_order_picks_the_same_subunit();
    test_ranked_rank_rules();
    test_ranked_declines_symmetric();
    test_within_molecule_symmetric_components();
    test_chain_order_invariant();
    test_lr_chain_refines_all_the_way();
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
