// Observable values at session build (issue #65).
//
// Session build used to settle obs_values by full-walking every
// observable, then build the incremental tracker's tables by computing
// the same per-molecule embedding counts a second time.  It now builds
// the tracker first and settles the tracked observables out of its
// tables — the contribs ARE the walk's per-molecule terms, and the
// per-complex pass flags ARE its per-complex verdicts.
//
// So the values a run starts from are now derived rather than walked
// for, and this pins both halves of that:
//
//   * they are the right numbers, against hand-computed counts over
//     `init_obs_seed_model.bngl`'s seed species;
//   * they are the SAME numbers a from-scratch walk of the identical
//     initial pool produces, observable by observable.
//
// The second check is the one that generalizes, and it has to be made
// across two sessions: `get_observable_values()` full-walks and
// overwrites obs_values, so asking one session for both readings would
// destroy the seeded one before it could be compared.
//
// Every assertion is exact — all values are integer counts over a seed
// state, so no trajectory statistics are involved.

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

double initial_value(const rulemonkey::Result& r, const std::string& name) {
  for (size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == name)
      return r.observable_data[i].front();
  std::fprintf(stderr, "FAIL: observable not found: %s\n", name.c_str());
  ++g_failures;
  return -1.0;
}

// The seed species are `A(b,s~0) 50`, `A(b!1,s~1).B(a!1) 30`,
// `A(b!1,s~0).A(b!1,s~0) 7`, `B(a) 20`, `C() 10`, so every count below
// is arithmetic on those, not a measurement:
//
//   Atot   A()             50 + 30 + 2*7          = 94   (trivial pattern)
//   Abound A(b!+)          30 + 2*7               = 44   (single-mol, constrained)
//   AB     A(b!1).B(a!1)   30                     = 30   (2-mol/1-bond)
//   AorB   A() B()         94 + (30 + 20)         = 144  (two patterns, two types)
//   ABsp   A(b!1).B(a!1)   30 complexes           = 30   (Species, tracked)
//   Csp    C()             10 complexes           = 10   (Species, untracked)
//   Apair  A()==2          the 7 A-A dimers       = 7    (Species, quantity form)
void test_initial_values_are_the_seed_species(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(9);
  auto r = sim.simulate(0.0, 5.0, 5);

  const struct {
    const char* name;
    double want;
  } expected[] = {
      {"Atot", 94}, {"Abound", 44}, {"AB", 30},   {"AorB", 144},
      {"ABsp", 30}, {"Csp", 10},    {"Apair", 7},
  };
  for (const auto& e : expected) {
    double const got = initial_value(r, e.name);
    check(got == e.want, std::string("t=0 value of ") + e.name + " should be " +
                             std::to_string(e.want) + " (got " + std::to_string(got) + ")");
  }
  sim.destroy_session();
}

// Whatever session build seeds must equal what a from-scratch walk of
// the same initial pool answers — for the tracked observables, which are
// now derived from the tracker's tables, and for the untracked ones,
// which are still walked for.  Two sessions: the walking one would
// otherwise overwrite the seeded values it is meant to check.
void test_seeded_values_match_a_full_walk(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator seeded(xml);
  seeded.initialize(9);
  auto r = seeded.simulate(0.0, 5.0, 5); // row 0 carries the seeded values

  rulemonkey::RuleMonkeySimulator walked(xml);
  walked.initialize(9);
  auto full = walked.get_observable_values(); // from-scratch, same initial pool

  check(full.size() == r.observable_names.size(),
        "the two sessions should report the same number of observables (" +
            std::to_string(full.size()) + " vs " + std::to_string(r.observable_names.size()) + ")");

  for (size_t i = 0; i < r.observable_names.size() && i < full.size(); ++i) {
    check(r.observable_data[i].front() == full[i],
          "seeded value of " + r.observable_names[i] + " should equal the full walk (seeded " +
              std::to_string(r.observable_data[i].front()) + ", walk " + std::to_string(full[i]) +
              ")");
  }
  seeded.destroy_session();
  walked.destroy_session();
}

// The model's `C() -> 0` rule is priced by `gateC() = if(ABsp>10, 5, 0)`,
// a rate law reading a tracked observable, so session build settles a
// rule against a value seeding produced.  At 5/s over ten molecules and
// five seconds, survival of even one C has probability ~1e-10.
void test_a_rate_law_reads_a_seeded_observable(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(9);
  auto r = sim.simulate(0.0, 5.0, 5);

  check(initial_value(r, "Csp") == 10, "C should start at its seeded population of 10");
  for (size_t i = 0; i < r.observable_names.size(); ++i) {
    if (r.observable_names[i] != "Csp")
      continue;
    check(r.observable_data[i].back() == 0,
          "the ABsp-gated decay should have consumed every C by t=5 (got " +
              std::to_string(r.observable_data[i].back()) + ")");
  }
  sim.destroy_session();
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <init_obs_seed_model.xml>\n", argv[0]);
    return 2;
  }
  const std::string xml = argv[1];

  try {
    test_initial_values_are_the_seed_species(xml);
    test_seeded_values_match_a_full_walk(xml);
    test_a_rate_law_reads_a_seeded_observable(xml);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: session build seeds every observable at the value a from-scratch "
                       "walk of the same pool answers\n");
  return 0;
}
