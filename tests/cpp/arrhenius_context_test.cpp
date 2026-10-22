// Which Arrhenius rules RuleMonkey runs, and which it refuses (issue #109).
//
// RuleMonkey runs an Arrhenius rule one of two ways.  Where the load-time
// expansion prices it exactly, it is expanded into conventional rules with
// fixed rates; anywhere else it runs as itself on the exact path, with dG
// counted at each firing as BioNetGen's network counts it per reaction.  The
// BNG2-oracle models arrhenius_context and arrhenius_context_bind pin the
// rates of the exact path; this test pins the line around it.  A rule
// neither path runs must be refused with an Error that names the rule and
// why, and a model either path runs must load with no Arrhenius diagnostic
// at all.
//
// Fixtures (tests/cpp/arrhenius_context/):
//   ep_disconnected   an energy pattern joined only by `.`, A().C(), whose
//                     count any bond can change
//   refused_shapes    a state change, two bonds at once, a homodimer, a ring
//                     closure, and exclude constraints, one rule each
//   scaffold_ring     a binding rule on the exact path whose only pair sits
//                     in one complex: it never fires, under -no-bscb too,
//                     since a rule with two reactant patterns joins two
//                     complexes, as BNG2 applies it
//   bivalent          a binding rule on the exact path, run again with its
//                     AddBond's two sites named the other way round: the
//                     run must be the same, bit for bit
//   no_energy_patterns  a rule on the exact path in a model with no energy
//                     patterns: it runs at exp(-Ea0)
//   ep_disconnected_apart  a pattern joined only by `.` naming a molecule
//                     no bond can bring near the rule's: it loads clean
//   ep_disconnected_ring  a pattern joined only by `.`, within reach of an
//                     unbinding written with `.`, which never parts a
//                     complex, and asking nothing of its sites: clean
//   ep_disconnected_site  the same with a pattern that requires A's x free,
//                     which the unbinding changes: refused
// Models that must load clean, passed after the fixture directory: the
// expansion's exact class and the exact path's oracle models.

#include "rulemonkey/simulator.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <regex>
#include <sstream>
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

// The Arrhenius diagnostics a model loads with.
std::vector<rulemonkey::UnsupportedFeature> arrhenius_features(const std::string& xml) {
  std::vector<rulemonkey::UnsupportedFeature> out;
  rulemonkey::RuleMonkeySimulator const sim(xml);
  for (const auto& f : sim.unsupported_features())
    if (f.element == "RateLaw@type=Arrhenius")
      out.push_back(f);
  return out;
}

struct Refusal {
  std::string rule;   // BNGL rule name
  std::string reason; // the phrase that must appear in its diagnostic
};

// `xml` must refuse exactly the rules in `want`, each with an Error whose text
// names the rule and carries its reason.
void expect_refusals(const std::string& label, const std::string& xml,
                     const std::vector<Refusal>& want) {
  std::vector<rulemonkey::UnsupportedFeature> got;
  try {
    got = arrhenius_features(xml);
  } catch (const std::exception& e) {
    check(false, label + ": load threw: " + e.what());
    return;
  }
  check(got.size() == want.size(), label + ": expected " + std::to_string(want.size()) +
                                       " Arrhenius refusal(s), got " + std::to_string(got.size()));
  for (const auto& w : want) {
    bool found = false;
    for (const auto& f : got) {
      if (!contains(f.feature, "rule '" + w.rule + "'"))
        continue;
      found = true;
      check(f.severity == rulemonkey::Severity::Error,
            label + ": rule " + w.rule + " must be refused with an Error, not a Warn");
      check(contains(f.feature, w.reason), label + ": rule " + w.rule +
                                               " refused without the reason \"" + w.reason +
                                               "\": " + f.feature);
    }
    check(found, label + ": rule " + w.rule + " was not refused");
  }
}

// The observables of one run.
std::vector<std::vector<double>> run_observables(const std::string& xml, double t_end,
                                                 std::uint64_t seed, bool bscb) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_block_same_complex_binding(bscb);
  return sim.run(rulemonkey::TimeSpec{0.0, t_end, 4}, seed).observable_data;
}

// A copy of `xml` with every AddBond's two sites swapped, written beside the
// build's temporary files; returns its path.
std::string with_addbond_swapped(const std::string& xml) {
  std::ifstream const in(xml);
  std::stringstream buf;
  buf << in.rdbuf();
  std::string const text = buf.str();
  std::regex const add_bond(R"re(<AddBond site1="([^"]*)" site2="([^"]*)")re");
  std::string const swapped =
      std::regex_replace(text, add_bond, R"re(<AddBond site1="$2" site2="$1")re");
  check(swapped != text, xml + " has no AddBond to swap");
  std::string const out = xml + ".swapped.xml";
  std::ofstream(out) << swapped;
  return out;
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <arrhenius_context dir> [clean.xml ...]\n", argv[0]);
    return 2;
  }
  std::string const dir = std::string(argv[1]) + "/";

  expect_refusals("disconnected pattern", dir + "ep_disconnected.xml",
                  {{"Unbind", "energy pattern A().C() is not connected"}});
  expect_refusals("refused shapes", dir + "refused_shapes.xml",
                  {{"Flip", "it changes a state or adds or deletes a molecule"},
                   {"Two", "it makes or breaks more than one bond"},
                   {"Homo", "its bond joins two molecules of one type"},
                   {"Ring", "it does not bind two reactant patterns to each other"},
                   {"Excl", "it carries include/exclude constraints"}});

  expect_refusals("disconnected pattern out of reach", dir + "ep_disconnected_apart.xml", {});
  expect_refusals("disconnected pattern, ring opening", dir + "ep_disconnected_ring.xml", {});
  expect_refusals("disconnected pattern on a bond site", dir + "ep_disconnected_site.xml",
                  {{"Open", "energy pattern A(x).C() is not connected"}});
  // A(b!+).X() never matches, so it adds nothing to the A-B unbinding.
  {
    double sum = 0;
    double sum2 = 0;
    int const n = 40;
    for (std::uint64_t seed = 1; seed <= static_cast<std::uint64_t>(n); ++seed) {
      auto const obs = run_observables(dir + "ep_disconnected_apart.xml", 1.0, seed, true);
      double const v = obs.empty() ? 0 : obs[0].back();
      sum += v;
      sum2 += v * v;
    }
    double const mean = sum / n;
    double const se = std::sqrt((sum2 / n - mean * mean) / (n - 1));
    double const want = 100.0 * std::exp(-std::exp(-0.9));
    check(std::abs(mean - want) < 6 * se, "disconnected pattern out of reach: A-B left " +
                                              std::to_string(mean) + " +- " + std::to_string(se) +
                                              ", want " + std::to_string(want));
  }

  // With no energy patterns, dG is 0: A(s~1).B unbinds at exp(-Ea0) and
  // A(s~0).B never.
  expect_refusals("no energy patterns", dir + "no_energy_patterns.xml", {});
  {
    double sum = 0;
    double sum2 = 0;
    int const n = 20;
    for (std::uint64_t seed = 1; seed <= static_cast<std::uint64_t>(n); ++seed) {
      auto const obs = run_observables(dir + "no_energy_patterns.xml", 1.0, seed, true);
      check(obs.size() == 2 && obs[1].back() == 1000.0,
            "no energy patterns: an A(s~0).B unbound, seed " + std::to_string(seed));
      double const v = obs.empty() ? 0 : obs[0].back();
      sum += v;
      sum2 += v * v;
    }
    double const mean = sum / n;
    double const se = std::sqrt((sum2 / n - mean * mean) / (n - 1));
    double const want = 1000.0 * std::exp(-std::exp(-0.5));
    check(std::abs(mean - want) < 6 * se, "no energy patterns: A(s~1).B left " +
                                              std::to_string(mean) + " +- " + std::to_string(se) +
                                              ", want " + std::to_string(want));
  }

  // The scaffold's A and B are in one complex, so the A-B bond never forms.
  for (std::uint64_t seed = 1; seed <= 20; ++seed) {
    for (bool const bscb : {true, false}) {
      auto const obs = run_observables(dir + "scaffold_ring.xml", 10.0, seed, bscb);
      check(!obs.empty() && obs[0].back() == 0.0,
            std::string("scaffold_ring: an A-B pair in one complex bound") +
                (bscb ? "" : " under -no-bscb") + ", seed " + std::to_string(seed));
    }
  }

  // The bivalent run is the same whichever order its AddBond names its sites.
  std::string const swapped = with_addbond_swapped(dir + "bivalent.xml");
  for (std::uint64_t seed = 1; seed <= 3; ++seed)
    check(run_observables(dir + "bivalent.xml", 5.0, seed, true) ==
              run_observables(swapped, 5.0, seed, true),
          "bivalent: swapping the AddBond's sites changed the run, seed " + std::to_string(seed));
  std::remove(swapped.c_str());

  for (int i = 2; i < argc; ++i)
    expect_refusals(std::string("clean ") + argv[i], argv[i], {});

  if (g_failures != 0) {
    std::fprintf(stderr, "arrhenius_context_test: %d check(s) failed\n", g_failures);
    return 1;
  }
  std::printf("arrhenius_context_test: OK\n");
  return 0;
}
