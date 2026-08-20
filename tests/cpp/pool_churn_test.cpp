// Pool bookkeeping under churn (issue #62).
//
// Three AgentPool/observable structures were rewritten to stop paying a
// per-population cost on every event or every sample.  Each swapped an
// O(N) walk for an O(1) side table, so each can now go quietly wrong
// instead of merely slowly:
//
//   * `type_mol_index_` removal is a swap-with-back keyed by a position
//     table.  A stale position drops a live molecule out of its type
//     index, or strands a dead one in it.
//   * `active_molecule_count()` is a maintained tally.  Drift moves the
//     population at which a `-gml` molecule limit stops the run.
//   * the Species full-walk dedupes complexes by stamping their member
//     molecules rather than by hashing complex ids, and short-circuits
//     the complex lookup for unbonded molecules.  A missed stamp counts
//     one complex once per member.
//
// Every assertion here is exact.  `pool_churn_model.bngl` conserves
// Xn + Yn and Zmol by construction, and its only molecule-count-changing
// rule is a zero-order synthesis, so the population a limit stops at is
// a fixed integer rather than a trajectory statistic.

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

const std::vector<double>& series(const rulemonkey::Result& r, const std::string& name) {
  for (size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == name)
      return r.observable_data[i];
  static const std::vector<double> empty;
  std::fprintf(stderr, "FAIL: observable not found: %s\n", name.c_str());
  ++g_failures;
  return empty;
}

std::string at(const std::string& name, size_t i, double v) {
  return name + "[" + std::to_string(i) + "]=" + std::to_string(v);
}

// X() <-> Y() is a molecule-type change: every fire deletes one molecule
// and adds another, which is exactly the traffic the type index sees.
// Xn + Yn can only stay at X0 if every live molecule is in exactly one
// type index and every dead one is in none.
void test_type_index_survives_churn(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(11);
  auto r = sim.simulate(0.0, 60.0, 60);

  const auto& xn = series(r, "Xn");
  const auto& yn = series(r, "Yn");
  const auto& zmol = series(r, "Zmol");

  bool saw_conversion = false;
  for (size_t i = 0; i < xn.size(); ++i) {
    check(xn[i] + yn[i] == 200.0, "X<->Y conserves the pool: " + at("Xn", i, xn[i]) + " + " +
                                      at("Yn", i, yn[i]) + " should be 200");
    check(zmol[i] == 100.0, "Z is never created or destroyed: " + at("Zmol", i, zmol[i]));
    if (yn[i] > 0.0)
      saw_conversion = true;
  }
  check(saw_conversion,
        "the X->Y rule should have fired at least once (test is otherwise vacuous)");

  // The same tally read through the public per-type accessor, which
  // walks the type index directly.
  int const live_x = sim.get_molecule_count("X");
  int const live_y = sim.get_molecule_count("Y");
  check(live_x + live_y == 200, "get_molecule_count over X and Y should total 200 (got " +
                                    std::to_string(live_x) + " + " + std::to_string(live_y) + ")");
  check(sim.get_molecule_count("Z") == 100, "get_molecule_count(Z) should be 100 (got " +
                                                std::to_string(sim.get_molecule_count("Z")) + ")");
  sim.destroy_session();
}

// `Zsp` counts complexes holding at least one Z; `Zmol` counts Z
// molecules; `Zdim` counts the two-Z complexes.  A dimer therefore costs
// Zmol two and Zsp one, so Zmol - Zsp is exactly the dimer count — an
// identity that breaks the moment the walk counts a complex once per
// member.  `Zpair` asks the same question through the quantity branch
// (`Z()==2`), which additionally requires the walk to have summed the
// match count across BOTH members of the complex, not just the one it
// entered through.
void test_species_walk_dedupes_complexes(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(23);
  auto r = sim.simulate(0.0, 60.0, 60);

  const auto& zmol = series(r, "Zmol");
  const auto& zsp = series(r, "Zsp");
  const auto& zdim = series(r, "Zdim");
  const auto& zpair = series(r, "Zpair");

  bool saw_dimer = false;
  for (size_t i = 0; i < zmol.size(); ++i) {
    check(zmol[i] - zsp[i] == zdim[i],
          "one complex per dimer, two molecules: " + at("Zmol", i, zmol[i]) + " - " +
              at("Zsp", i, zsp[i]) + " should equal " + at("Zdim", i, zdim[i]));
    check(zpair[i] == zdim[i], "the `Z()==2` quantity form should agree with the bond form: " +
                                   at("Zpair", i, zpair[i]) + " vs " + at("Zdim", i, zdim[i]));
    if (zdim[i] > 0.0)
      saw_dimer = true;
  }
  check(saw_dimer, "dimers should have formed (test is otherwise vacuous)");
  sim.destroy_session();
}

// The pool starts at 300 molecules (200 X + 100 Z) and only `0 -> W()`
// changes that total, one molecule at a time.  The engine stops the run
// on the first event that puts the count ABOVE the limit, so a limit of
// 350 must halt at 351 molecules — 51 W's — and not one more.
//
// The breaching event is not rolled back, so the trajectory's last row
// and the pool left behind must agree: `Wn` is an incrementally tracked
// observable and the trailing sample rows are emitted after the loop
// exits, which is exactly where a row can end up describing the pool as
// it was one event earlier.
void test_molecule_limit_stops_at_the_cap(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_molecule_limit(350);
  sim.initialize(5);
  auto r = sim.simulate(0.0, 100.0, 20);

  const auto& wn = series(r, "Wn");
  check(wn.back() == 51.0,
        "a 350-molecule limit over a 300-molecule pool should stop at 51 synthesized W (got " +
            std::to_string(wn.back()) + ")");
  check(sim.get_molecule_count("W") == 51, "the session pool should hold 51 W (got " +
                                               std::to_string(sim.get_molecule_count("W")) + ")");
  int const total = sim.get_molecule_count("X") + sim.get_molecule_count("Y") +
                    sim.get_molecule_count("Z") + sim.get_molecule_count("W");
  check(total == 351,
        "the run should stop at exactly limit+1 molecules (got " + std::to_string(total) + ")");
  sim.destroy_session();

  // Without the limit the same horizon runs far past that population —
  // otherwise the check above would pass on a limit that never engaged.
  rulemonkey::RuleMonkeySimulator unbounded(xml);
  unbounded.initialize(5);
  auto ur = unbounded.simulate(0.0, 100.0, 20);
  check(series(ur, "Wn").back() > 51.0, "the unlimited run should synthesize well past 51 W (got " +
                                            std::to_string(series(ur, "Wn").back()) + ")");
  unbounded.destroy_session();
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <pool_churn_model.xml>\n", argv[0]);
    return 2;
  }
  const std::string xml = argv[1];

  try {
    test_type_index_survives_churn(xml);
    test_species_walk_dedupes_complexes(xml);
    test_molecule_limit_stops_at_the_cap(xml);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: type indices, the live-molecule tally, and the Species complex walk "
                       "all stay exact under add/delete churn\n");
  return 0;
}
