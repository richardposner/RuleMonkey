// Load-time diagnostics for the two MM(kcat,Km) constructs where RM cannot
// reproduce BNG2 (issue #45).
//
// BNG2.pl is the reference RM is written against, so both of these are
// divergences from it. They warn rather than refuse because the constructs
// are idiomatic BNGL and refusing would put a large share of real MM models
// out of reach. What the warning buys is that the divergence is named at
// load rather than discovered by diffing trajectories against BNG2.
//
// (a) A reactant pattern that can match more than one species. BNG2's
//     network expansion emits one MM reaction per matching species PAIR,
//     each evaluating the law on that pair's own counts, while RM applies
//     one law to the summed match counts. Measured: 2.00x for a two-species
//     substrate in saturation, 1.81x for a two-species enzyme with the
//     enzyme in excess, the law being nonlinear in both counts. Matching
//     BNG2 would need a live species-canonical-form to match-count map on
//     each slot, maintained per event, which is the species bookkeeping a
//     network-free engine exists to avoid. RM checks both slots, BNG2 having
//     run its own check on the substrate alone. BNG2 warns at rule-read time
//     via checkSpeciesGraph(..., IsSpecies => 1); RM recomputes the
//     predicate, and this test pins that the two agree on that axis —
//     the predicate was validated against BNG2 2.9.3 on complete patterns,
//     a missing component, an unspecified state, `!+` and `!?` bonds, a
//     stateless molecule type and a complete two-molecule complex, agreeing
//     on all seven.
//
// (b) A symmetry_factor that cannot be attributed. It belongs to the
//     reactant pattern the rule transforms (issue #37); when the rule
//     transforms BOTH, the scalar is a product of the two patterns' own
//     factors and the XML is one number with no way to split it. RM applies
//     it to the substrate, which reproduces BNG2 for the canonical shape and
//     anywhere the law is linear, and runs up to 2x fast against BNG2 in
//     saturation if the symmetry was the enzyme slot's.
//
// Fixtures:
//   mm_diagnostics_model.xml  R1 trips (a) on the substrate axis, R3 trips
//                             (a) on the enzyme axis, R2 trips (b). Each
//                             rule trips exactly one: R2's own patterns are
//                             deliberately complete species.
//   mm_symmetry_model.xml     negative control: four MM rules, including the
//                             symmetric-substrate ones of #37 and the
//                             enzyme-slot one whose factor IS attributable.
//                             None of them may warn.
//
// Both are Warn severity, so a model carrying them still runs; the CLI's
// exit code is unaffected.

#include "rulemonkey/simulator.hpp"

#include <cstdio>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

bool contains(const std::string& haystack, const std::string& needle) {
  return haystack.find(needle) != std::string::npos;
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <mm_diagnostics_model.xml> <mm_symmetry_model.xml>\n", argv[0]);
    return 2;
  }
  std::string const diag_xml = argv[1];
  std::string const clean_xml = argv[2];

  try {
    rulemonkey::RuleMonkeySimulator const sim(diag_xml);
    const auto& uf = sim.unsupported_features();

    int substrate_axis = 0;
    int enzyme_axis = 0;
    int attribution = 0;
    for (const auto& f : uf) {
      bool const is_substrate = contains(f.feature, "whose substrate pattern can match");
      bool const is_enzyme = contains(f.feature, "whose enzyme pattern can match");
      bool const is_attribution = contains(f.feature, "transforms both of its reactant patterns");
      if (is_substrate || is_enzyme || is_attribution)
        check(f.severity == rulemonkey::Severity::Warn,
              "MM diagnostics must be Warn: the constructs are idiomatic BNGL and refusing "
              "them would put a large share of real MM models out of reach");
      // Matched by BNGL rule name: BNG2 does not preserve rule order in the
      // XML ids (R3 lands on RR2 here).
      if (is_substrate) {
        ++substrate_axis;
        check(contains(f.feature, "(R1)"),
              "substrate-axis diagnostic must name rule R1, got: " + f.feature.substr(0, 90));
      }
      if (is_enzyme) {
        ++enzyme_axis;
        check(contains(f.feature, "(R3)"),
              "enzyme-axis diagnostic must name rule R3, got: " + f.feature.substr(0, 90));
      }
      if (is_attribution) {
        ++attribution;
        check(contains(f.feature, "(R2)"),
              "attribution diagnostic must name rule R2, got: " + f.feature.substr(0, 90));
      }
      // Every rule in the fixture trips exactly one of the three.
      check(static_cast<int>(is_substrate) + static_cast<int>(is_enzyme) +
                    static_cast<int>(is_attribution) <=
                1,
            "one entry tripped more than one diagnostic");
    }
    check(substrate_axis == 1,
          "expected exactly 1 substrate-axis diagnostic, got " + std::to_string(substrate_axis));
    check(enzyme_axis == 1,
          "expected exactly 1 enzyme-axis diagnostic, got " + std::to_string(enzyme_axis));
    check(attribution == 1,
          "expected exactly 1 attribution diagnostic, got " + std::to_string(attribution));

    // Nothing here is Error severity: both constructs run.
    for (const auto& f : uf)
      check(f.severity != rulemonkey::Severity::Error,
            "no Error-severity feature expected on this fixture, got: " + f.feature.substr(0, 80));

    std::printf("mm_diagnostics_model: %d substrate-axis, %d enzyme-axis, %d attribution, "
                "%zu features\n",
                substrate_axis, enzyme_axis, attribution, uf.size());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: loading %s threw: %s\n", diag_xml.c_str(), e.what());
    return 1;
  }

  // Negative control: four MM rules RM does reproduce BNG2 on.
  try {
    rulemonkey::RuleMonkeySimulator const sim(clean_xml);
    const auto& uf = sim.unsupported_features();
    for (const auto& f : uf) {
      check(!contains(f.feature, "more than one species"),
            "mm_symmetry_model tripped the multi-species diagnostic: " + f.feature.substr(0, 120));
      check(!contains(f.feature, "transforms both of its reactant patterns"),
            "mm_symmetry_model tripped the attribution diagnostic — its RD arm transforms "
            "only the enzyme slot, so the factor IS attributable: " +
                f.feature.substr(0, 120));
    }
    std::printf("mm_symmetry_model: %zu features, none of them an MM diagnostic\n", uf.size());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: loading %s threw: %s\n", clean_xml.c_str(), e.what());
    return 1;
  }

  if (g_failures) {
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
  }
  std::printf("mm_diagnostics_test: all checks passed\n");
  return 0;
}
