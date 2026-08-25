// A rate law over a Species observable must not full-walk the pool on
// every event (issue #79).
//
// A Species observable that a rate law reads was excluded from the
// incremental tracker and served instead by
// compute_rate_dependent_observables — a from-scratch walk of every
// complex holding a molecule of the pattern's seed type, run after
// EVERY SSA event rather than once per sample.  So a rule as small as
//
//     Species Rtot R()
//     rate()  = kt*Rtot
//     R() -> P()   rate()
//
// paid O(pool) per event and the whole run went quadratic in the seed
// population: 70 / 177 / 511 us per event at 10k / 30k / 100k R, against
// a flat ~3 us for the same chemistry whose rate law happens not to read
// an observable.  The exclusion was deliberate — a tracked Species obs
// only settles when its dirty complexes are flushed, and the flush fired
// at sample time, far too late for a rate law read on the next event.
// The fix flushes those obs (and only those) after every event, which is
// O(complexes this event dirtied) rather than O(pool).
//
// Two arms, because the defect has two halves:
//
// ARM 1 — cost, and the freshness that makes the cheap path legal.
//   dynamic_rate_species_obs_model    Species   Rtot R()
//   dynamic_rate_molecules_obs_model  Molecules Rtot R()
//
//   Every R is a free monomer, so the two observables hold the same
//   number at every instant of the run; the Molecules one was already
//   maintained per-event by the well-tested per-mid delta path.  Same
//   arithmetic on the same RNG stream, and the only difference is which
//   path keeps the value the rate law reads.  That makes the Molecules
//   twin both the cost floor and an exact trajectory oracle:
//
//     1. same seed -> identical trajectories.  This is the freshness
//        assertion.  A rate law reading a value one event stale draws
//        against a different propensity and the two series separate —
//        which is exactly what happens with the per-event flush removed
//        (checked: t=10 lands on 22974 instead of 23048).
//     2. same seed -> comparable wall time.  This is the defect.  The
//        ratio is ~1.6 with the fix and ~60 without, so the bound below
//        is loose by a factor of ~6 in either direction and does not
//        depend on how fast the machine running it is.
//
// ARM 2 — the bookkeeping the per-event flush leans on.
//   dynamic_rate_species_cx_model declares one Species pattern twice:
//   `Chains` is read by a rate law and so takes the new per-event flush,
//   `ChainsCtl` is not and so keeps the once-per-sample cadence.  They
//   share the per-molecule contribution tables and the pool's
//   dead-complex side channel — which is drained, not copied — so a
//   per-event flush that swallowed a notification the per-sample obs
//   still needed would show up here as the two disagreeing.  A(s,t)
//   polymerises into chains, so bond ops move non-endpoint molecules
//   between complexes rather than only capping at dimers.
//
//   The independent oracle is get_observable_values(), which full-walks
//   every observable from scratch.  Nothing during a run ever re-derives
//   a tracked value that way, so any delta the tracker got wrong at any
//   event of the run is still there at the end to be caught.

#include "rulemonkey/simulator.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
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

int idx_of(const std::vector<std::string>& names, const std::string& name) {
  for (size_t i = 0; i < names.size(); ++i)
    if (names[i] == name)
      return static_cast<int>(i);
  return -1;
}

// Wall time of run() alone — the constructor's XML parse is identical
// work for both arms and would only dilute the ratio.
double timed_run(const std::string& xml, std::uint64_t seed, const rulemonkey::TimeSpec& ts,
                 rulemonkey::Result& out) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  auto const t0 = std::chrono::steady_clock::now();
  out = sim.run(ts, seed);
  auto const t1 = std::chrono::steady_clock::now();
  return std::chrono::duration<double>(t1 - t0).count();
}

const std::vector<double>& series(const rulemonkey::Result& r, const std::string& name) {
  static const std::vector<double> empty;
  int const i = idx_of(r.observable_names, name);
  if (i < 0) {
    check(false, "observable '" + name + "' missing from result");
    return empty;
  }
  return r.observable_data[static_cast<size_t>(i)];
}

const std::vector<double>& fn_series(const rulemonkey::Result& r, const std::string& name) {
  static const std::vector<double> empty;
  int const i = idx_of(r.function_names, name);
  if (i < 0) {
    check(false, "function '" + name + "' missing from result");
    return empty;
  }
  return r.function_data[static_cast<size_t>(i)];
}

constexpr std::uint64_t kSeed = 7;

// ---- ARM 1 -------------------------------------------------------------

void test_species_obs_matches_molecules_twin(const std::string& species_xml,
                                             const std::string& molecules_xml) {
  const rulemonkey::TimeSpec ts{0.0, 10.0, 10};
  rulemonkey::Result r_species, r_molecules;
  double const t_species = timed_run(species_xml, kSeed, ts, r_species);
  double const t_molecules = timed_run(molecules_xml, kSeed, ts, r_molecules);

  const auto& s_obs = series(r_species, "Rtot");
  const auto& m_obs = series(r_molecules, "Rtot");
  const auto& s_fn = fn_series(r_species, "conv_rate");
  const auto& m_fn = fn_series(r_molecules, "conv_rate");

  check(s_obs.size() == m_obs.size() && !s_obs.empty(),
        "both arms should produce the same number of sample rows");

  // Exact, not tolerance-based: identical arithmetic on an identical RNG
  // stream.  A stale read shifts a propensity and the streams separate.
  bool obs_same = s_obs.size() == m_obs.size();
  bool fn_same = s_fn.size() == m_fn.size();
  for (size_t i = 0; i < s_obs.size() && i < m_obs.size(); ++i)
    if (s_obs[i] != m_obs[i])
      obs_same = false;
  for (size_t i = 0; i < s_fn.size() && i < m_fn.size(); ++i)
    if (s_fn[i] != m_fn[i])
      fn_same = false;

  check(obs_same, "Species-obs and Molecules-obs arms should give identical Rtot trajectories "
                  "(a rate law reading a stale Species obs separates them)");
  check(fn_same, "Species-obs and Molecules-obs arms should give identical conv_rate values");

  if (!obs_same) {
    for (size_t i = 0; i < s_obs.size() && i < m_obs.size(); ++i)
      std::fprintf(stderr, "  row %2zu: species=%.17g molecules=%.17g\n", i, s_obs[i], m_obs[i]);
  }

  // The defect.  Floor is ~1.6x; without the fix this arm is ~60x.
  constexpr double kMaxRatio = 10.0;
  double const ratio = t_species / (t_molecules > 0 ? t_molecules : 1e-9);
  check(ratio < kMaxRatio,
        "a rate law over a Species obs should cost about what the same rate law over a "
        "Molecules obs costs; ratio was " +
            std::to_string(ratio) + " (species " + std::to_string(t_species) + "s, molecules " +
            std::to_string(t_molecules) + "s)");
  std::fprintf(stderr, "arm 1: species %.4fs, molecules %.4fs, ratio %.2f\n", t_species,
               t_molecules, ratio);
}

// ---- ARM 2 -------------------------------------------------------------

void test_species_cx_bookkeeping(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  auto const r = sim.run({0.0, 60.0, 12}, kSeed);

  const auto& chains = series(r, "Chains");
  const auto& ctl = series(r, "ChainsCtl");
  const auto& rate = fn_series(r, "chain_rate");
  double const kt = sim.get_parameter("kt");

  check(!chains.empty() && chains.size() == ctl.size(),
        "Chains and ChainsCtl should have the same number of sample rows");

  // The two obs declare the same pattern; only the flush cadence differs.
  bool agree = true;
  for (size_t i = 0; i < chains.size() && i < ctl.size(); ++i) {
    if (chains[i] != ctl[i]) {
      agree = false;
      std::fprintf(stderr, "  row %2zu: Chains=%.17g ChainsCtl=%.17g\n", i, chains[i], ctl[i]);
    }
  }
  check(agree, "the rate-read Species obs (per-event flush) and its unread twin "
               "(per-sample flush) must hold the same value at every sample");

  // The rate law's own column must be the observable it reads, times kt.
  bool rate_ok = rate.size() == chains.size();
  for (size_t i = 0; i < rate.size() && i < chains.size(); ++i)
    if (std::fabs(rate[i] - (kt * chains[i])) > 1e-12 * (1.0 + std::fabs(rate[i])))
      rate_ok = false;
  check(rate_ok, "chain_rate column should equal kt*Chains at every sample");

  // Independent oracle: a from-scratch walk of every observable, taken
  // once at t_end on a live session driven to the same instant with the
  // same seed.  Nothing during a run ever re-derives a tracked value
  // that way, so a delta the tracker got wrong at any event of the run
  // is still there to be caught — and the call is made only after the
  // last step, because it re-bases obs_values and would heal any drift
  // it did not first report.
  //
  // Every observable is checked, not just the Species pair: A_bound and
  // Q_tot ride the long-settled Molecules delta path, so if they agree
  // and the Species pair does not, the tracker is at fault rather than
  // the trajectory.
  rulemonkey::RuleMonkeySimulator live(xml);
  live.initialize(kSeed);
  live.step_to(60.0);
  auto const truth = live.get_observable_values();
  const auto& names = live.observable_names();
  for (const char* nm : {"Chains", "ChainsCtl", "A_bound", "Q_tot"}) {
    int const i = idx_of(names, nm);
    check(i >= 0, std::string("observable '") + nm + "' should exist");
    if (i < 0 || chains.empty())
      continue;
    double const tracked = series(r, nm).back();
    check(truth[static_cast<size_t>(i)] == tracked,
          std::string("incrementally tracked '") + nm + "' (" + std::to_string(tracked) +
              ") should equal a from-scratch walk (" +
              std::to_string(truth[static_cast<size_t>(i)]) + ") at t_end");
  }

  // Sanity: the model has to actually exercise complexes and the rule
  // the observable prices, or the assertions above are vacuous.
  check(chains.back() > 0, "the model should have formed chains by t_end");
  check(series(r, "Q_tot").back() > 0, "the rate-priced rule should have fired");
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 4) {
    std::fprintf(stderr, "usage: %s <species_obs.xml> <molecules_obs.xml> <species_cx.xml>\n",
                 argv[0]);
    return 2;
  }
  try {
    test_species_obs_matches_molecules_twin(argv[1], argv[2]);
    test_species_cx_bookkeeping(argv[3]);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: unexpected exception: %s\n", e.what());
    return 1;
  }
  if (g_failures != 0) {
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "all checks passed\n");
  return 0;
}
