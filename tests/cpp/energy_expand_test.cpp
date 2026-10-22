// Direct unit tests for the eBNGL Sekar energy-rule expansion
// (cpp/rulemonkey/energy_expand.{hpp,cpp}).
//
// The oracle is the reference NFsim (RuleWorld/nfsim, commit c4f1bb2), whose
// eBNGL loader prints the pre-computed rate constants for each expanded
// variant.  For the test-suite model `v40` (A(b)+B(a) <-> A(b!1).B(a!1) with
// energy pattern A(b!1,x!2).B(a!1).X(a!2), phi=0.5, RT=1.0, Ea0=2.0, g_val=5.0)
// NFsim prints exactly:
//
//   Expanding _R1 with 1 context condition(s), 2 variant(s) per direction:
//     v0: ΔG=0  k_fwd=0.135335  k_rev=0.135335  context=[A.x=free]
//     v1: ΔG=5  k_fwd=0.011109  k_rev=1.64872   context=[A.x=bound]
//
// These hardcoded values are a non-circular oracle: they come from a
// separately-built binary, not from our own exp() calls.  At RT = 1 NFsim's
// rates for a reversible binding rule are BioNetGen's, so they also check
// arrhenius_rate, which prices each rule from the variant's ΔG (issue #107).
// A second case (the cooperative-scaffold model from nfsim
// validate/basicModels/v40) pins the multi-pattern always/conditional
// classification and ΔG summing.  A third pins which rule takes phi and which
// 1-phi, against the rates BioNetGen's network gives the rules of
// tests/bng_oracle/models/arrhenius_formula.bngl.

#include "energy_expand.hpp"

#include <cmath>
#include <cstdio>
#include <string>
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

void check_close(double got, double want, double tol, const std::string& msg) {
  if (std::fabs(got - want) > tol) {
    std::fprintf(stderr, "FAIL: %s (got %.6g, want %.6g)\n", msg.c_str(), got, want);
    ++g_failures;
  }
}

// Build a bound component with an optional intra-pattern partner recorded.
EpComponent bound_comp(const std::string& name) {
  EpComponent c;
  c.name = name;
  c.is_bound = true;
  return c;
}

EpComponent free_comp(const std::string& name) {
  EpComponent c;
  c.name = name;
  c.is_free = true;
  return c;
}

// A component the pattern lists but asks nothing of (`!?`).
EpComponent any_comp(const std::string& name) {
  EpComponent c;
  c.name = name;
  return c;
}

// ---- Case 1: v40 test-suite model — single conditional context pattern ----
void test_v40_conditional_context() {
  // Energy pattern A(b!1,x!2).B(a!1).X(a!2) = g_val (5.0).
  EnergyPatternInfo ep;
  ep.id = "EP1";
  ep.energy_value = 5.0;
  ep.energy_expr = "g_val"; // symbolic source, for the re-resolvable rate
  ep.molecules = {
      EpMolecule{"A", {bound_comp("b"), bound_comp("x")}},
      EpMolecule{"B", {bound_comp("a")}},
      EpMolecule{"X", {bound_comp("a")}},
  };
  ep.bonds = {
      {0, 0, 1, 0}, // A.b — B.a  (reaction center)
      {0, 1, 2, 0}, // A.x — X.a  (extra context)
  };

  EnergyFunction efn;
  efn.add_pattern(ep);

  // Reaction center A.b — B.a, phi=0.5, Ea0=2.0.
  auto variants = efn.expand_binding("A", "b", "B", "a");

  check(variants.size() == 2, "v40: expected 2 context variants");
  if (variants.size() != 2)
    return;

  // Locate variants by their single context constraint on A.x.
  const EnergyExpandedVariant* free_v = nullptr;
  const EnergyExpandedVariant* bound_v = nullptr;
  for (const auto& v : variants) {
    check(v.constraints.size() == 1, "v40: variant must carry one context constraint");
    if (v.constraints.empty())
      continue;
    check(v.constraints[0].reactant_idx == 0, "v40: A.x constraint is on reactant 0 (A)");
    check(v.constraints[0].comp_name == "x", "v40: context constraint is on component x");
    if (v.constraints[0].must_be_bound)
      bound_v = &v;
    else
      free_v = &v;
  }
  check(free_v != nullptr && bound_v != nullptr, "v40: need one free and one bound variant");
  if (free_v == nullptr || bound_v == nullptr)
    return;

  // v0 (A.x free): ΔG = 0, symbolic ΔG is "0".
  check_close(free_v->delta_g, 0.0, 1e-9, "v40 v0: ΔG");
  check(free_v->delta_g_expr == "0", "v40 v0: symbolic ΔG is \"0\"");
  // The binding rule takes phi on the binding ΔG; its reverse half takes
  // 1-phi on the unbinding ΔG.
  check_close(arrhenius_rate(2.0, 0.5, free_v->delta_g, false), 0.135335, 1e-5, "v40 v0: k_fwd");
  check_close(arrhenius_rate(2.0, 0.5, -free_v->delta_g, true), 0.135335, 1e-5, "v40 v0: k_rev");

  // v1 (A.x bound): ΔG = 5, symbolic ΔG references the source parameter.
  check_close(bound_v->delta_g, 5.0, 1e-9, "v40 v1: ΔG");
  check(bound_v->delta_g_expr == "(g_val)", "v40 v1: symbolic ΔG is \"(g_val)\"");
  check_close(arrhenius_rate(2.0, 0.5, bound_v->delta_g, false), 0.011109, 1e-5, "v40 v1: k_fwd");
  check_close(arrhenius_rate(2.0, 0.5, -bound_v->delta_g, true), 1.64872, 1e-4, "v40 v1: k_rev");
}

// ---- Case 2: cooperative scaffold — always + conditional patterns ---------
void test_cooperative_scaffold() {
  // From nfsim validate/basicModels/v40.bngl:
  //   ep1: S(A!1).A(s!1)             Gf_SA  = -47.5   (always, for S.A rule)
  //   ep2: S(B!1).B(s!1)             Gf_SB  = -47.5   (irrelevant to S.A rule)
  //   ep3: S(A!1,B!2).A(s!1).B(s!2)  Gf_SAB =   1.8   (conditional: S.B bound)
  // Rule S(A)+A(s): reaction center S.A — A.s.
  const double gf_sa = -47.5;
  const double gf_sb = -47.5;
  const double gf_sab = 1.8;

  EnergyPatternInfo ep1;
  ep1.id = "EP1";
  ep1.energy_value = gf_sa;
  ep1.molecules = {EpMolecule{"S", {bound_comp("A")}}, EpMolecule{"A", {bound_comp("s")}}};
  ep1.bonds = {{0, 0, 1, 0}}; // S.A — A.s

  EnergyPatternInfo ep2;
  ep2.id = "EP2";
  ep2.energy_value = gf_sb;
  ep2.molecules = {EpMolecule{"S", {bound_comp("B")}}, EpMolecule{"B", {bound_comp("s")}}};
  ep2.bonds = {{0, 0, 1, 0}}; // S.B — B.s

  EnergyPatternInfo ep3;
  ep3.id = "EP3";
  ep3.energy_value = gf_sab;
  ep3.molecules = {EpMolecule{"S", {bound_comp("A"), bound_comp("B")}},
                   EpMolecule{"A", {bound_comp("s")}}, EpMolecule{"B", {bound_comp("s")}}};
  ep3.bonds = {{0, 0, 1, 0}, {0, 1, 2, 0}}; // S.A — A.s, S.B — B.s

  EnergyFunction efn;
  efn.add_pattern(ep1);
  efn.add_pattern(ep2);
  efn.add_pattern(ep3);

  auto variants = efn.expand_binding("S", "A", "A", "s");

  check(variants.size() == 2, "coop: expected 2 context variants (S.B free / bound)");
  if (variants.size() != 2)
    return;

  const EnergyExpandedVariant* free_v = nullptr;
  const EnergyExpandedVariant* bound_v = nullptr;
  for (const auto& v : variants) {
    check(v.constraints.size() == 1, "coop: one context constraint (S.B)");
    if (v.constraints.empty())
      continue;
    check(v.constraints[0].reactant_idx == 0, "coop: S.B constraint on reactant 0 (S)");
    check(v.constraints[0].comp_name == "B", "coop: context constraint on component B");
    if (v.constraints[0].must_be_bound)
      bound_v = &v;
    else
      free_v = &v;
  }
  check(free_v != nullptr && bound_v != nullptr, "coop: need free and bound variants");
  if (free_v == nullptr || bound_v == nullptr)
    return;

  // S.B free: only the always pattern ep1 contributes → ΔG = Gf_SA.
  check_close(free_v->delta_g, gf_sa, 1e-6, "coop v0: ΔG = Gf_SA");
  // S.B bound: ep1 + ep3 → ΔG = Gf_SA + Gf_SAB.
  check_close(bound_v->delta_g, gf_sa + gf_sab, 1e-6, "coop v1: ΔG = Gf_SA + Gf_SAB");
}

// ---- Case 3: no overlapping energy pattern → ΔG = 0, single variant -------
void test_no_relevant_pattern() {
  EnergyPatternInfo ep; // pattern with an unrelated bond
  ep.id = "EP1";
  ep.energy_value = 3.0;
  ep.molecules = {EpMolecule{"C", {bound_comp("y")}}, EpMolecule{"D", {bound_comp("z")}}};
  ep.bonds = {{0, 0, 1, 0}};

  EnergyFunction efn;
  efn.add_pattern(ep);
  auto variants = efn.expand_binding("A", "b", "B", "a");
  check(variants.size() == 1, "no-overlap: single variant");
  if (variants.empty())
    return;
  check(variants[0].constraints.empty(), "no-overlap: no context constraints");
  check_close(variants[0].delta_g, 0.0, 1e-9, "no-overlap: ΔG = 0");
  check_close(arrhenius_rate(1.5, 0.5, variants[0].delta_g, false), std::exp(-1.5), 1e-12,
              "no-overlap: k = exp(-Ea0), no RT");
}

// ---- Case 4: phi or 1-phi, and no RT --------------------------------------
// The rates BioNetGen's network gives the rules of
// tests/bng_oracle/models/arrhenius_formula.bngl (phi = 0.3), each with dG the
// free-energy change of the reaction the rule writes.
void test_bng2_phi_formula() {
  double const phi = 0.3;
  double const ln10 = std::log(10.0);
  // Break: Va(b!1).Vb(a!1) -> Va(b) + Vb(a), unidirectional; the binding dG
  // is -ln(10), so the unbinding dG is +ln(10), and the rule takes phi.
  check_close(arrhenius_rate(6.216979751083923, phi, ln10, false), 1e-3, 1e-12,
              "unidirectional bond-breaking rule: phi on its own dG");
  // Revu: Qa(b!1).Qb(a!1) <-> Qa(b) + Qb(a), binding dG ln(10).  The rule as
  // written unbinds (dG -ln(10)) and takes phi; its reverse half binds (dG
  // ln(10)) and takes 1-phi.
  check_close(arrhenius_rate(9.901115899874398, phi, -ln10, false), 1e-4, 1e-13,
              "reversible rule written as unbinding: phi on the unbinding");
  check_close(arrhenius_rate(9.901115899874398, phi, ln10, true), 1e-5, 1e-14,
              "its reverse half: 1-phi on the binding");
  // Rev: Ra(b) + Rb(a) <-> Ra(b!1).Rb(a!1), binding dG ln(100).
  check_close(arrhenius_rate(10.131374409173802, phi, 2.0 * ln10, false), 1e-5, 1e-14,
              "reversible binding rule: phi on the binding");
  check_close(arrhenius_rate(10.131374409173802, phi, -2.0 * ln10, true), 1e-3, 1e-12,
              "its reverse half: 1-phi on the unbinding");
}

// ---- Case 5: which patterns the expansion prices exactly (issue #109) -----
// inexact_reason against the reaction center A.b-B.a, one pattern at a time.
// The issue's nine reproducers are pinned end to end by
// arrhenius_context_test; these are the branches they do not reach, and the
// shapes that must pass.
std::string reason_for(const EnergyPatternInfo& ep, const SiteGraph& sites = {},
                       const std::string& t1 = "A", const std::string& s1 = "b",
                       const std::string& t2 = "B", const std::string& s2 = "a") {
  EnergyFunction efn;
  efn.add_pattern(ep);
  return efn.inexact_reason(t1, s1, t2, s2, sites);
}

EnergyPatternInfo pattern(const std::string& text, std::vector<EpMolecule> mols,
                          std::vector<EnergyPatternInfo::Bond> bonds) {
  EnergyPatternInfo ep;
  ep.id = "EP1";
  ep.text = text;
  ep.energy_value = 1.0;
  ep.molecules = std::move(mols);
  ep.bonds = std::move(bonds);
  return ep;
}

void test_inexact_reason() {
  auto expect = [](const std::string& got, const std::string& want, const std::string& label) {
    if (want.empty())
      check(got.empty(), label + ": must expand, refused: " + got);
    else
      check(got.find(want) != std::string::npos,
            label + ": want \"" + want + "\", got \"" + got + "\"");
  };

  // Shapes that must pass.
  expect(reason_for(pattern("A(b!1).B(a!1)", {{"A", {bound_comp("b")}}, {"B", {bound_comp("a")}}},
                            {{0, 0, 1, 0}})),
         "", "plain pair");
  expect(reason_for(pattern("B(a!1).A(b!1)", {{"B", {bound_comp("a")}}, {"A", {bound_comp("b")}}},
                            {{0, 0, 1, 0}})),
         "", "plain pair, written B first");
  expect(reason_for(pattern("A(b!1,c!?).B(a!1)",
                            {{"A", {bound_comp("b"), any_comp("c")}}, {"B", {bound_comp("a")}}},
                            {{0, 0, 1, 0}})),
         "", "a site asked nothing of");
  EnergyPatternInfo const ctx = pattern("A(b!1,c!2).B(a!1).C(a!2)",
                                        {{"A", {bound_comp("b"), bound_comp("c")}},
                                         {"B", {bound_comp("a")}},
                                         {"C", {bound_comp("a")}}},
                                        {{0, 0, 1, 0}, {0, 1, 2, 0}});
  expect(reason_for(ctx), "", "context bond, A.c never bound");
  SiteGraph only_c;
  only_c.partners[{"A", "c"}] = {{"C", "a"}};
  expect(reason_for(ctx, only_c), "", "context bond, A.c binds only C.a");
  expect(reason_for(pattern("A(b!1).B(a!1,c!2).C(a!2)",
                            {{"A", {bound_comp("b")}},
                             {"B", {bound_comp("a"), bound_comp("c")}},
                             {"C", {bound_comp("a")}}},
                            {{0, 0, 1, 0}, {1, 1, 2, 0}})),
         "", "context bond on B");
  expect(reason_for(pattern("A(b!1).C(a!1)", {{"A", {bound_comp("b")}}, {"C", {bound_comp("a")}}},
                            {{0, 0, 1, 0}})),
         "", "A.b bonded to another partner: count unchanged");
  expect(reason_for(pattern("A(c)", {{"A", {free_comp("c")}}}, {})), "",
         "a free site off the reaction center: count unchanged");

  // Branches the issue's reproducers do not reach.
  expect(reason_for(pattern("A(b).B(a)", {{"A", {any_comp("b")}}, {"B", {any_comp("a")}}}, {})),
         "energy pattern A(b).B(a) is not connected", "disconnected");
  expect(reason_for(pattern("A(b)", {{"A", {free_comp("b")}}}, {})),
         "energy pattern A(b) constrains A.b without naming the A.b-B.a bond",
         "reaction site required free");
  expect(reason_for(pattern("B(a!+)", {{"B", {bound_comp("a")}}}, {})),
         "energy pattern B(a!+) constrains B.a without naming the A.b-B.a bond",
         "reaction site bound to anything");
  expect(reason_for(pattern("A(b!1).B(a!1,a!2).A(b!2)",
                            {{"A", {bound_comp("b")}},
                             {"B", {bound_comp("a"), bound_comp("a")}},
                             {"A", {bound_comp("b")}}},
                            {{0, 0, 1, 0}, {1, 1, 2, 0}})),
         "names the A.b-B.a bond more than once", "center bond twice");
  expect(reason_for(pattern(
             "A(b!1,c!2).B(a!1,d!2)",
             {{"A", {bound_comp("b"), bound_comp("c")}}, {"B", {bound_comp("a"), bound_comp("d")}}},
             {{0, 0, 1, 0}, {0, 1, 1, 1}})),
         "closes a ring with the A.b-B.a bond", "ring");
  expect(reason_for(pattern("A(b!1,c!2).B(a!1).A(x!2)",
                            {{"A", {bound_comp("b"), bound_comp("c")}},
                             {"B", {bound_comp("a")}},
                             {"A", {bound_comp("x")}}},
                            {{0, 0, 1, 0}, {0, 1, 2, 0}})),
         "names a second A, a molecule type of the reaction center", "third molecule of type A");
  SiteGraph repeated;
  repeated.repeated[{"A", "c"}] = 2;
  expect(reason_for(ctx, repeated), "requires A.c bound, but A declares more than one c",
         "repeated context site");
  SiteGraph never;
  never.partners[{"A", "c"}] = {{"D", "a"}};
  expect(reason_for(ctx, never), "requires A.c bound to C.a, but A.c binds D.a",
         "context pairing that never occurs");

  // A lone `!+` on a center molecule is the expansion's own condition.
  EnergyPatternInfo const plus = pattern(
      "A(b!1,c!+).B(a!1)", {{"A", {bound_comp("b"), bound_comp("c")}}, {"B", {bound_comp("a")}}},
      {{0, 0, 1, 0}});
  expect(reason_for(plus), "", "lone !+ on a center molecule");
  expect(reason_for(plus, repeated), "requires A.c bound, but A declares more than one c",
         "lone !+ on a repeated site");
  expect(reason_for(pattern(
             "A(b!1,c!+,d!+).B(a!1)",
             {{"A", {bound_comp("b"), bound_comp("c"), bound_comp("d")}}, {"B", {bound_comp("a")}}},
             {{0, 0, 1, 0}})),
         "requires more than one site besides the A.b-B.a bond bound", "two !+");
  expect(reason_for(pattern("A(b!1,c!2,d!+).B(a!1).C(a!2)",
                            {{"A", {bound_comp("b"), bound_comp("c"), bound_comp("d")}},
                             {"B", {bound_comp("a")}},
                             {"C", {bound_comp("a")}}},
                            {{0, 0, 1, 0}, {0, 1, 2, 0}})),
         "requires more than one site besides the A.b-B.a bond bound", "!+ and a context bond");

  // A `!?` counts once only where it can land one way.
  EnergyPatternInfo const loose =
      pattern("A(b!1).B(a!1,x!?)",
              {{"A", {bound_comp("b")}}, {"B", {bound_comp("a"), any_comp("x")}}}, {{0, 0, 1, 0}});
  SiteGraph two_x;
  two_x.repeated[{"B", "x"}] = 2;
  expect(reason_for(loose, two_x), "can place its B.x (!?) 2 ways", "!? on a repeated site");
  EnergyPatternInfo const tri =
      pattern("L(r!1,r!?).R(l!1)",
              {{"L", {bound_comp("r"), any_comp("r")}}, {"R", {bound_comp("l")}}}, {{0, 0, 1, 0}});
  SiteGraph two_r;
  two_r.repeated[{"L", "r"}] = 2;
  expect(reason_for(tri, two_r, "L", "r", "R", "l"), "", "!? on the one other r of L(r,r)");
  SiteGraph three_r;
  three_r.repeated[{"L", "r"}] = 3;
  expect(reason_for(tri, three_r, "L", "r", "R", "l"), "can place its L.r (!?) 2 ways",
         "!? on L(r,r,r)");
  // Two interchangeable `!?` on L(r,r,r)'s two other r: 2!/0! matches over
  // 2! automorphisms, once.
  EnergyPatternInfo const tri2 =
      pattern("L(r!1,r!?,r!?).R(l!1)",
              {{"L", {bound_comp("r"), any_comp("r"), any_comp("r")}}, {"R", {bound_comp("l")}}},
              {{0, 0, 1, 0}});
  expect(reason_for(tri2, three_r, "L", "r", "R", "l"), "", "two !? on L(r,r,r)");
}

} // namespace

int main() {
  test_v40_conditional_context();
  test_cooperative_scaffold();
  test_no_relevant_pattern();
  test_bng2_phi_formula();
  test_inexact_reason();

  if (g_failures != 0) {
    std::fprintf(stderr, "energy_expand_test: %d check(s) failed\n", g_failures);
    return 1;
  }
  std::printf("energy_expand_test: all checks passed\n");
  return 0;
}
