// Unit tests for the exact Arrhenius free-energy machinery
// (cpp/rulemonkey/energy_match.{hpp,cpp}, issue #109).
//
// The oracle is BioNetGen's own definition: an energy pattern weighs
// matches / automorphisms in a species (EnergyPattern.pm), the
// automorphisms are the pattern's matches on itself (SpeciesGraph.pm), and a
// reaction's dG is the sum of G times the change in weight from reactants to
// products.  Each expected dG below is that sum worked by hand for the
// species before and after, and the cases are the issue's: the dG the issue
// quotes for each reproducer's reaction is the one asserted here.  The
// matcher is run on a hand-built graph, so these pin the counting alone,
// with no engine and no expansion.

#include "energy_match.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

using namespace rulemonkey;

namespace {

int g_failures = 0;

void check(bool cond, const std::string& msg) {
  if (!cond) {
    std::fprintf(stderr, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

void check_close(double got, double want, const std::string& msg) {
  if (std::fabs(got - want) > 1e-12) {
    std::fprintf(stderr, "FAIL: %s (got %.12g, want %.12g)\n", msg.c_str(), got, want);
    ++g_failures;
  }
}

// Molecule types and component names, by number.
// Plain ints by design: they index the toy graph and the patterns directly.
enum Type { A, B, C, D, L, R, S, kTypes };       // NOLINT(cppcoreguidelines-use-enum-class)
enum Name { b, a, c, d, r, l, s, x, m, kNames }; // NOLINT(cppcoreguidelines-use-enum-class)

// A graph of molecules and components, read through the matcher's view.
struct Toy {
  std::vector<int> type;
  std::vector<int> comp_mol, comp_name, comp_state, comp_partner;

  // A molecule with the given (name, state) components; returns its
  // components' ids in order.
  std::vector<int> mol(int t, const std::vector<std::pair<int, int>>& comps) {
    int const m = static_cast<int>(type.size());
    type.push_back(t);
    std::vector<int> ids;
    for (auto [name, state] : comps) {
      ids.push_back(static_cast<int>(comp_mol.size()));
      comp_mol.push_back(m);
      comp_name.push_back(name);
      comp_state.push_back(state);
      comp_partner.push_back(-1);
    }
    return ids;
  }
  void bond(int c1, int c2) {
    comp_partner[c1] = c2;
    comp_partner[c2] = c1;
  }

  int type_of(int m) const { return type[m]; }
  int mol_of(int comp) const { return comp_mol[comp]; }
  int name_of(int comp) const { return comp_name[comp]; }
  int state_of(int comp) const { return comp_state[comp]; }
  int partner_of(int comp) const { return comp_partner[comp]; }
  void comps_named(int m, int name, std::vector<int>& out) const {
    for (int k = 0; k < static_cast<int>(comp_mol.size()); ++k)
      if (comp_mol[k] == m && comp_name[k] == name)
        out.push_back(k);
  }
};

// Pattern components.
EmComponent to(int name, int mol, int comp, int state = -1) {
  EmComponent e;
  e.name = name;
  e.state = state;
  e.bond = EmComponent::Bond::To;
  e.to_mol = mol;
  e.to_comp = comp;
  return e;
}
EmComponent plus(int name) {
  EmComponent e;
  e.name = name;
  e.bond = EmComponent::Bond::Bound;
  return e;
}
EmComponent any(int name, int state = -1) {
  EmComponent e;
  e.name = name;
  e.state = state;
  return e;
}
EmComponent none(int name) {
  EmComponent e;
  e.name = name;
  e.bond = EmComponent::Bond::Free;
  return e;
}

EmPattern pat(std::vector<EmMolecule> mols, double g) {
  EmPattern p;
  p.mols = std::move(mols);
  p.energy = g;
  p.automorphisms = em_automorphisms(p);
  return p;
}

// The model's patterns, all relevant.
ArrheniusModel model_of(std::vector<EmPattern> ps) {
  ArrheniusModel am;
  am.patterns = std::move(ps);
  // Slot names for the rooted bounds: each type's components, as the toy
  // molecules declare them.
  am.slot_names.resize(kTypes);
  am.slot_names[A] = {b, c, d, s};
  am.slot_names[B] = {a, c, d, s};
  am.slot_names[C] = {a, b, d};
  am.slot_names[D] = {a, c};
  am.slot_names[L] = {r, r, r};
  am.slot_names[R] = {l};
  am.slot_names[S] = {l};
  return am;
}

std::vector<int> all(const ArrheniusModel& am) {
  std::vector<int> v;
  v.reserve(am.patterns.size());
  for (int i = 0; i < static_cast<int>(am.patterns.size()); ++i)
    v.push_back(i);
  return v;
}

// A(b!1).B(a!1) — the pair every case below starts from, worth G_AB = -2.
EmPattern pair_ab() { return pat({{A, {to(b, 1, 0)}}, {B, {to(a, 0, 0)}}}, -2.0); }
EmPattern pair_lr() { return pat({{L, {to(r, 1, 0)}}, {R, {to(l, 0, 0)}}}, -2.0); }

void test_automorphisms() {
  check_close(pair_ab().automorphisms, 1, "A(b!1).B(a!1) has one automorphism");
  check_close(pat({{L, {plus(r), plus(r)}}}, 0).automorphisms, 2, "L(r!+,r!+): its two r swap");
  check_close(pat({{L, {to(r, 1, 0), to(r, 2, 0)}}, {R, {to(l, 0, 0)}}, {R, {to(l, 0, 1)}}}, 0)
                  .automorphisms,
              2, "L(r!1,r!2).R(l!1).R(l!2): the two R swap");
  check_close(pat({{L, {to(r, 1, 0), any(r), any(r)}}, {R, {to(l, 0, 0)}}}, 0).automorphisms, 2,
              "L(r!1,r!?,r!?).R(l!1): the two r!? swap");
  check_close(pat({{L, {to(r, 1, 0), plus(r)}}, {R, {to(l, 0, 0)}}}, 0).automorphisms, 1,
              "L(r!1,r!+).R(l!1): nothing swaps");
  // BioNetGen counts a component's bonds, a `!+` or `!?` as one, so a `!+`
  // lands on a `!?` and L(r!+,r!?) has two automorphisms.  A `~?` needs some
  // state, so it does not land on a component with none: L(r~?!+,r!+) has one.
  check_close(pat({{L, {plus(r), any(r)}}}, 0).automorphisms, 2, "L(r!+,r!?): the two r swap");
  EmComponent q = plus(r);
  q.state = kEmAnyState;
  check_close(pat({{L, {q, plus(r)}}}, 0).automorphisms, 1, "L(r~?!+,r!+): nothing swaps");
}

// Breaking the A-B bond of each case's complexes: the dG the issue quotes.
void test_unbinding_cases() {
  // Case 1: A(b!1).B(a!1,s~1) G_S1 = -1.
  {
    ArrheniusModel const am =
        model_of({pair_ab(), pat({{A, {to(b, 1, 0)}}, {B, {to(a, 0, 0), any(s, 1)}}}, -1.0)});
    Toy g;
    auto ab1 = g.mol(A, {{b, -1}});
    auto bb1 = g.mol(B, {{a, -1}, {s, 1}});
    g.bond(ab1[0], bb1[0]);
    auto ab0 = g.mol(A, {{b, -1}});
    auto bb0 = g.mol(B, {{a, -1}, {s, 0}});
    g.bond(ab0[0], bb0[0]);
    check_close(em_delta_g_bond(am, all(am), g, ab1[0], bb1[0], false), 3, "case 1, s~1");
    check_close(em_delta_g_bond(am, all(am), g, ab0[0], bb0[0], false), 2, "case 1, s~0");
  }
  // Case 2: A(b!1,c).B(a!1) G_free = -1.
  {
    ArrheniusModel const am =
        model_of({pair_ab(), pat({{A, {to(b, 1, 0), none(c)}}, {B, {to(a, 0, 0)}}}, -1.0)});
    Toy g;
    auto a1 = g.mol(A, {{b, -1}, {c, -1}});
    auto b1 = g.mol(B, {{a, -1}});
    g.bond(a1[0], b1[0]);
    auto a2 = g.mol(A, {{b, -1}, {c, -1}});
    auto b2 = g.mol(B, {{a, -1}});
    auto c2 = g.mol(C, {{a, -1}});
    g.bond(a2[0], b2[0]);
    g.bond(a2[1], c2[0]);
    check_close(em_delta_g_bond(am, all(am), g, a1[0], b1[0], false), 3, "case 2, c free");
    check_close(em_delta_g_bond(am, all(am), g, a2[0], b2[0], false), 2, "case 2, c bound");
  }
  // Case 3: L(r!1,r!2).R(l!1).R(l!2) with L(r!+,r!+) G_full = -1.
  {
    ArrheniusModel const am = model_of({pair_lr(), pat({{L, {plus(r), plus(r)}}}, -1.0)});
    Toy g;
    auto lm = g.mol(L, {{r, -1}, {r, -1}, {r, -1}});
    auto r1 = g.mol(R, {{l, -1}});
    auto r2 = g.mol(R, {{l, -1}});
    g.bond(lm[0], r1[0]);
    g.bond(lm[1], r2[0]);
    check_close(em_delta_g_bond(am, all(am), g, lm[0], r1[0], false), 3, "case 3, doubly bound L");
  }
  // Case 4: L(r!1,r!+).R(l!1) G_full = -1, on a singly bound L(r,r).
  {
    ArrheniusModel const am =
        model_of({pair_lr(), pat({{L, {to(r, 1, 0), plus(r)}}, {R, {to(l, 0, 0)}}}, -1.0)});
    Toy g;
    auto lm = g.mol(L, {{r, -1}, {r, -1}});
    auto r1 = g.mol(R, {{l, -1}});
    g.bond(lm[0], r1[0]);
    check_close(em_delta_g_bond(am, all(am), g, lm[0], r1[0], false), 2,
                "case 4, the other r free");
  }
  // Case 5: L(r!1,r!2).R(l!1).S(l!2) G_RS = -1, on L.R.S and L.R.R.
  {
    ArrheniusModel const am = model_of(
        {pair_lr(),
         pat({{L, {to(r, 1, 0), to(r, 2, 0)}}, {R, {to(l, 0, 0)}}, {S, {to(l, 0, 1)}}}, -1.0)});
    Toy g;
    auto l1 = g.mol(L, {{r, -1}, {r, -1}});
    auto rr = g.mol(R, {{l, -1}});
    auto ss = g.mol(S, {{l, -1}});
    g.bond(l1[0], rr[0]);
    g.bond(l1[1], ss[0]);
    auto l2 = g.mol(L, {{r, -1}, {r, -1}});
    auto ra = g.mol(R, {{l, -1}});
    auto rb = g.mol(R, {{l, -1}});
    g.bond(l2[0], ra[0]);
    g.bond(l2[1], rb[0]);
    check_close(em_delta_g_bond(am, all(am), g, l1[0], rr[0], false), 3, "case 5, L.R.S");
    check_close(em_delta_g_bond(am, all(am), g, l2[0], ra[0], false), 2, "case 5, L.R.R");
  }
  // Case 8: A(b!1,c!2,d!3).B(a!1).C(a!2).D(a!3) G_CD = -1.
  {
    ArrheniusModel const am =
        model_of({pair_ab(), pat({{A, {to(b, 1, 0), to(c, 2, 0), to(d, 3, 0)}},
                                  {B, {to(a, 0, 0)}},
                                  {C, {to(a, 0, 1)}},
                                  {D, {to(a, 0, 2)}}},
                                 -1.0)});
    Toy g;
    auto a1 = g.mol(A, {{b, -1}, {c, -1}, {d, -1}});
    auto b1 = g.mol(B, {{a, -1}});
    auto c1 = g.mol(C, {{a, -1}});
    auto d1 = g.mol(D, {{a, -1}});
    g.bond(a1[0], b1[0]);
    g.bond(a1[1], c1[0]);
    g.bond(a1[2], d1[0]);
    auto a2 = g.mol(A, {{b, -1}, {c, -1}, {d, -1}});
    auto b2 = g.mol(B, {{a, -1}});
    auto c2 = g.mol(C, {{a, -1}});
    g.bond(a2[0], b2[0]);
    g.bond(a2[1], c2[0]);
    check_close(em_delta_g_bond(am, all(am), g, a1[0], b1[0], false), 3, "case 8, A.B.C.D");
    check_close(em_delta_g_bond(am, all(am), g, a2[0], b2[0], false), 2, "case 8, A.B.C");
  }
  // A `!?` counted C(m, k) times: L(r!1,r!?).R(l!1) on L(r,r,r) twice per
  // bond, L(r!1,r!?,r!?).R(l!1) once.
  {
    ArrheniusModel const am =
        model_of({pair_lr(), pat({{L, {to(r, 1, 0), any(r)}}, {R, {to(l, 0, 0)}}}, -1.0),
                  pat({{L, {to(r, 1, 0), any(r), any(r)}}, {R, {to(l, 0, 0)}}}, -0.25)});
    Toy g;
    auto lm = g.mol(L, {{r, -1}, {r, -1}, {r, -1}});
    auto r1 = g.mol(R, {{l, -1}});
    g.bond(lm[0], r1[0]);
    check_close(em_delta_g_bond(am, all(am), g, lm[0], r1[0], false), 2 + (2 * 1) + 0.25,
                "!? on L(r,r,r)");
  }
}

// Making a bond between two complexes: the direct count against its cut
// into one term per end and a product across the bond.
void test_binding_cut() {
  // A(b!1,c!2).B(a!1,d!3).C(a!2).D(a!3) G = -1: context on both sides.
  EmPattern const both = pat({{A, {to(b, 1, 0), to(c, 2, 0)}},
                              {B, {to(a, 0, 0), to(d, 3, 0)}},
                              {C, {to(a, 0, 1)}},
                              {D, {to(a, 1, 1)}}},
                             -1.0);
  // A(b!+) G = 0.5: a term on A's end alone.
  EmPattern const a_bound = pat({{A, {plus(b)}}}, 0.5);
  ArrheniusModel am = model_of({pair_ab(), both, a_bound});
  Toy g;
  auto am1 = g.mol(A, {{b, -1}, {c, -1}});
  auto cm = g.mol(C, {{a, -1}});
  g.bond(am1[1], cm[0]);
  auto bm = g.mol(B, {{a, -1}, {d, -1}});
  auto dm = g.mol(D, {{a, -1}});
  g.bond(bm[1], dm[0]);
  double const direct = em_delta_g_bond(am, all(am), g, am1[0], bm[0], true);
  check_close(direct, -2 - 1 + 0.5, "binding, direct count");

  // The same from the cut: the end terms plus, per cut bond, the product of
  // the two sides' rooted matches over the pattern's automorphisms.
  double cut_sum = em_delta_g_end(am, all(am), g, am1[0]) + em_delta_g_end(am, all(am), g, bm[0]);
  for (int pi = 0; pi < 2; ++pi) {
    const EmPattern& p = am.patterns[pi];
    EmPattern s0, s1;
    int m0 = -1, c0 = -1, m1 = -1, c1 = -1;
    bool const ok = em_cut(p, 0, 0, s0, m0, c0, s1, m1, c1);
    check(ok, "cut of pattern " + std::to_string(pi));
    if (!ok)
      continue;
    long const left = em_count_rooted(s0, g, m0, c0, am1[0]);
    long const right = em_count_rooted(s1, g, m1, c1, bm[0]);
    cut_sum += p.energy * static_cast<double>(left * right) / p.automorphisms;
    check(em_rooted_bound(s1, m1, c1, am) >= static_cast<double>(right),
          "side 1's bound covers its matches, pattern " + std::to_string(pi));
  }
  check_close(cut_sum, direct, "binding, cut count equals direct count");

  // A bond on a cycle cannot be cut.
  EmPattern const ring = pat({{A, {to(b, 1, 0), to(c, 1, 1)}}, {B, {to(a, 0, 0), to(c, 0, 1)}}}, 0);
  EmPattern s0, s1;
  int m0 = -1, c0 = -1, m1 = -1, c1 = -1;
  check(!em_cut(ring, 0, 0, s0, m0, c0, s1, m1, c1), "a bond on a cycle is not cut");

  check(em_connected(both), "a bonded pattern is connected");
  check(!em_connected(pat({{A, {any(b)}}, {B, {any(a)}}}, 0)), "A().B() is not connected");
}

// The matcher's own rules, on shapes the cases above do not reach.
void test_matching_rules() {
  // A pattern with a cycle, A(b!1,c!2).B(a!1,d!2): both of A's bonds must
  // land on the same B, which a match built bond by bond from A checks only
  // once every molecule is placed.
  EmPattern const ring = pat({{A, {to(b, 1, 0), to(c, 1, 1)}}, {B, {to(a, 0, 0), to(d, 0, 1)}}}, 0);
  Toy g;
  auto a1 = g.mol(A, {{b, -1}, {c, -1}});
  auto b1 = g.mol(B, {{a, -1}, {d, -1}});
  g.bond(a1[0], b1[0]);
  g.bond(a1[1], b1[1]);
  // A chain B3-A2-B2-A3: every site the ring names is bound, but A2's c
  // holds another B than its b does.
  auto a2 = g.mol(A, {{b, -1}, {c, -1}});
  auto b2 = g.mol(B, {{a, -1}, {d, -1}});
  auto b3 = g.mol(B, {{a, -1}, {d, -1}});
  auto a3 = g.mol(A, {{b, -1}, {c, -1}});
  g.bond(a2[0], b2[0]);
  g.bond(a2[1], b3[1]);
  g.bond(a3[1], b2[1]);
  check(em_count_rooted(ring, g, 0, 0, a1[0]) == 1, "a ring matches the ring");
  check(em_count_rooted(ring, g, 0, 0, a2[0]) == 0, "a ring does not match a chain");

  // Two pattern molecules never land on one molecule:
  // L(r!1).R(l!1,m!2).L(r!2) needs two L, so an R holding both r of one L
  // does not match it.
  EmPattern const two_l =
      pat({{L, {to(r, 1, 0)}}, {R, {to(l, 0, 0), to(m, 2, 0)}}, {L, {to(r, 1, 1)}}}, 0);
  Toy h;
  auto l1 = h.mol(L, {{r, -1}, {r, -1}});
  auto r1 = h.mol(R, {{l, -1}, {m, -1}});
  h.bond(l1[0], r1[0]);
  h.bond(l1[1], r1[1]);
  auto l2 = h.mol(L, {{r, -1}, {r, -1}});
  auto l3 = h.mol(L, {{r, -1}, {r, -1}});
  auto r2 = h.mol(R, {{l, -1}, {m, -1}});
  h.bond(l2[0], r2[0]);
  h.bond(l3[0], r2[1]);
  check(em_count_rooted(two_l, h, 1, 0, r1[0]) == 0, "two pattern L never share one L");
  check(em_count_rooted(two_l, h, 1, 0, r2[0]) == 1, "two pattern L on two L");
}

} // namespace

int main() {
  test_automorphisms();
  test_unbinding_cases();
  test_binding_cut();
  test_matching_rules();
  if (g_failures != 0) {
    std::fprintf(stderr, "energy_match_test: %d check(s) failed\n", g_failures);
    return 1;
  }
  std::printf("energy_match_test: all checks passed\n");
  return 0;
}
