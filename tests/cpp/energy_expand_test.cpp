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
// separately-built binary, not from our own exp() calls.  A second case
// (the cooperative-scaffold model from nfsim validate/basicModels/v40)
// pins the multi-pattern always/conditional classification and ΔG summing.

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

  EnergyFunction efn(/*rt=*/1.0);
  efn.add_pattern(ep);

  // Reaction center A.b — B.a, phi=0.5, Ea0=2.0.
  auto variants = efn.expand_binding("A", "b", "B", "a", /*ea0=*/2.0, /*phi=*/0.5);

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
  check_close(free_v->fwd_rate, 0.135335, 1e-5, "v40 v0: k_fwd");
  check_close(free_v->rev_rate, 0.135335, 1e-5, "v40 v0: k_rev");

  // v1 (A.x bound): ΔG = 5, symbolic ΔG references the source parameter.
  check_close(bound_v->delta_g, 5.0, 1e-9, "v40 v1: ΔG");
  check(bound_v->delta_g_expr == "(g_val)", "v40 v1: symbolic ΔG is \"(g_val)\"");
  check_close(bound_v->fwd_rate, 0.011109, 1e-5, "v40 v1: k_fwd");
  check_close(bound_v->rev_rate, 1.64872, 1e-4, "v40 v1: k_rev");
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
  const double rt = 2.5774863;
  const double phi = 0.5;
  const double ea0 = -17.8;

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

  EnergyFunction efn(rt);
  efn.add_pattern(ep1);
  efn.add_pattern(ep2);
  efn.add_pattern(ep3);

  auto variants = efn.expand_binding("S", "A", "A", "s", ea0, phi);

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

  // Rates follow k_fwd = exp(-(Ea0 + phi·ΔG)/RT); recompute independently.
  auto kf = [&](double dg) { return std::exp(-(ea0 + (phi * dg)) / rt); };
  auto kr = [&](double dg) { return std::exp(-(ea0 + ((phi - 1.0) * dg)) / rt); };
  check_close(free_v->fwd_rate, kf(gf_sa), 1e-6 * kf(gf_sa), "coop v0: k_fwd");
  check_close(free_v->rev_rate, kr(gf_sa), 1e-6 * kr(gf_sa), "coop v0: k_rev");
  check_close(bound_v->fwd_rate, kf(gf_sa + gf_sab), 1e-6 * kf(gf_sa + gf_sab), "coop v1: k_fwd");
  check_close(bound_v->rev_rate, kr(gf_sa + gf_sab), 1e-6 * kr(gf_sa + gf_sab), "coop v1: k_rev");
}

// ---- Case 3: no overlapping energy pattern → ΔG = 0, single variant -------
void test_no_relevant_pattern() {
  EnergyPatternInfo ep; // pattern with an unrelated bond
  ep.id = "EP1";
  ep.energy_value = 3.0;
  ep.molecules = {EpMolecule{"C", {bound_comp("y")}}, EpMolecule{"D", {bound_comp("z")}}};
  ep.bonds = {{0, 0, 1, 0}};

  EnergyFunction efn(2.0);
  efn.add_pattern(ep);
  auto variants = efn.expand_binding("A", "b", "B", "a", /*ea0=*/1.5, /*phi=*/0.5);
  check(variants.size() == 1, "no-overlap: single variant");
  if (variants.empty())
    return;
  check(variants[0].constraints.empty(), "no-overlap: no context constraints");
  check_close(variants[0].delta_g, 0.0, 1e-9, "no-overlap: ΔG = 0");
  check_close(variants[0].fwd_rate, std::exp(-1.5 / 2.0), 1e-9, "no-overlap: k_fwd = exp(-Ea0/RT)");
}

} // namespace

int main() {
  test_v40_conditional_context();
  test_cooperative_scaffold();
  test_no_relevant_pattern();

  if (g_failures != 0) {
    std::fprintf(stderr, "energy_expand_test: %d check(s) failed\n", g_failures);
    return 1;
  }
  std::printf("energy_expand_test: all checks passed\n");
  return 0;
}
