// Regression test for issue #24: rules with three or more reactant
// patterns must be reported, not silently skipped.
//
// The engine carries exactly two reactant slots.  Every consumer of
// `reactant_pattern_starts` treats slot B as the whole tail
// `[starts[1], molecules.size())`, so a rule written `A + A + A -> P`
// collapses its second and third patterns into one bond-free slot-B
// pattern.  `count_multi_mol_fast_generic` can only satisfy the
// unassigned molecules of such a pattern from inside the seed's own
// complex, so three free monomers score zero, `b_total` stays zero, and
// the rule's propensity is identically zero.  It never fires.
//
// The failure is invisible from the trajectory alone: mass is conserved
// and every other rule in the model behaves, so the run looks valid
// unless it is compared against another engine (NFsim fires the same
// rule ~126 times over the fixture's horizon).  This test pins down that
// RM refuses the model up front instead.
//
// The fixture (`trimolecular_model.bngl/.xml`) is the issue's reproducer:
// `r: A(s) + A(s) + A(s) -> P() k DeleteMolecules` with A(s) seeded at 400.
//
// NOTE: this test guards the *diagnostic*, which is a stopgap.  When the
// engine grows real n-ary reactant support the Tier-0 error goes away and
// this test must be replaced by a trajectory-parity check against NFsim.

#include "rulemonkey/simulator.hpp"

#include <cstdio>
#include <string>

namespace {

int g_failures = 0;

void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::fprintf(stdout, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stdout, "Usage: trimolecular_test <trimolecular_xml> <bimolecular_xml>\n");
    return 2;
  }
  const std::string tri_xml = argv[1];
  const std::string bi_xml = argv[2];

  // (1) A 3-reactant-pattern rule must raise an Error-severity feature,
  //     naming the rule and its pattern count.
  {
    const rulemonkey::RuleMonkeySimulator sim(tri_xml);
    const auto& uf = sim.unsupported_features();

    const rulemonkey::UnsupportedFeature* hit = nullptr;
    for (const auto& f : uf) {
      if (f.element == "ListOfReactantPatterns") {
        hit = &f;
        break;
      }
    }

    check(hit != nullptr, "3-reactant-pattern rule must report a ListOfReactantPatterns feature");
    if (hit != nullptr) {
      check(hit->severity == rulemonkey::Severity::Error,
            "3-reactant-pattern rule must be Error severity, not Warn — a Warn would "
            "still let the model run to a silently wrong trajectory");
      check(hit->feature.find("'r'") != std::string::npos,
            "diagnostic must name the offending rule so it is findable in a large model");
      check(hit->feature.find("3 reactant patterns") != std::string::npos,
            "diagnostic must state the pattern count it saw");
    }
  }

  // (2) The check must not fire on an ordinary bimolecular model — the
  //     cutoff is at three patterns, and two must stay clean.
  {
    const rulemonkey::RuleMonkeySimulator sim(bi_xml);
    for (const auto& f : sim.unsupported_features()) {
      check(f.element != "ListOfReactantPatterns",
            "bimolecular rule must not trip the >=3-reactant-pattern check");
    }
  }

  if (g_failures == 0)
    std::fprintf(stdout, "trimolecular_test: OK\n");
  return g_failures == 0 ? 0 : 1;
}
