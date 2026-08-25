// A rate law over a MULTI-PATTERN Species observable must not full-walk
// the pool on every event (issue #81).
//
// #79 / #80 removed the per-event full walk for a rate law reading a
// `Species` observable, by letting such an observable onto the
// incremental tracker and flushing its dirty complexes after each event.
// One shape was left behind: `init_incremental_observables` admitted a
// `Species` observable only when it declared exactly one pattern, so
//
//     Species Rtot  R(x~U), R(x~P), R(), S()
//     conv_rate() = kt*Rtot
//     R() -> P()  conv_rate()
//
// still fell through to compute_rate_dependent_observables and paid the
// same O(pool) walk per event that #79 measured on the single-pattern
// case.  The restriction was on `Species` alone — `Molecules` has taken
// any number of patterns all along.
//
// What made it more than a one-line gate is the semantics.  The full
// walk in `evaluate_observable` re-walks the pool once per pattern and
// adds one for every complex that pattern passes on its own, so a
// complex two of an observable's patterns match counts twice, and each
// pattern's quantity relation is applied to that pattern's own
// complex-wide count.  A single tally per complex — which is what the
// tracker kept — can represent neither.  So the per-complex tables are
// now per pattern, and these two arms are what pins that against the
// walk they have to agree with.
//
// ARM 1 — cost, and the freshness that makes the cheap path legal.
//   multi_pattern_species_obs_model    Species   Rtot R(x~U), R(x~P), R(), S()
//   multi_pattern_molecules_obs_model  Molecules Rtot R(x~U), R(x~P), R(), S()
//
//   Every R and S is a free monomer, so each pattern's species count is
//   that pattern's molecule count and the two observables hold the same
//   number at every instant — overlapping patterns included, since both
//   forms sum across patterns rather than merging them.  The Molecules
//   twin was already maintained per-event by the well-tested per-mid
//   delta path, so it is both the cost floor and an exact trajectory
//   oracle:
//
//     1. same seed -> identical trajectories.  This is the freshness
//        assertion: a rate law reading a value one event stale draws
//        against a different propensity and the two series separate.
//     2. same seed -> comparable wall time.  This is the defect.  The
//        ratio is ~1.7 with the fix; without it the Species arm does not
//        finish inside two minutes against the twin's 0.07 s.
//
// ARM 2 — the semantics, and the bookkeeping the per-event flush leans on.
//   multi_pattern_species_cx_model declares two overlapping-pattern
//   Species observables over a polymerising pool, both read by the rate
//   law and so both on the per-event flush:
//
//     Mixed   A()>=3, B()==1        different relations, different seed
//                                   types; both match a chain of >=3 A's
//                                   carrying exactly one B
//     Chains  A(s!+), A(s!+,t!+)    same seed type, one pattern a strict
//                                   subset of the other
//
//   The seeded trimer `A(s!1,t,b!3).A(s!2,t!1,b).A(s,t!2,b).B(a!3)`
//   makes t=0 checkable by hand: every copy contributes 2 to each, so
//   Mixed(0) = 2*T_tot + B_tot and Chains(0) = 2*T_tot.  One merged
//   tally per complex would give T_tot + B_tot and T_tot.  That row is
//   the tracker's own seeded value — `refresh_observables_for_sample`
//   leaves tracked observables alone — so it pins the seeding walk, and
//   `get_observable_values()` on a fresh session pins the full walk it
//   has to match.
//
//   Each has a twin declaring the same patterns and read by nothing, and
//   they are unread in two different ways.  `MixedWalk`'s patterns are
//   structurally unconstrained, so an unread observable of that shape is
//   not tracked at all and the engine full-walks it at every sample —
//   which makes it a from-scratch oracle for `Mixed` on every row rather
//   than only at t_end.  `ChainsCtl`'s patterns carry bond constraints,
//   so it IS tracked and keeps the once-per-sample flush cadence; it
//   shares the per-molecule contribution tables and the pool's
//   dead-complex side channel with `Chains`, and that channel is drained
//   rather than copied, so a per-event flush that swallowed a
//   notification the per-sample twin still needed shows up as the two
//   disagreeing.
//
//   The independent oracle at t_end is get_observable_values(), which
//   full-walks every observable from scratch.  Nothing during a run ever
//   re-derives a tracked value that way, so any delta the tracker got
//   wrong at any event of the run is still there at the end to be
//   caught.
//
// ARM 3 — the degenerate pattern the same classification has to survive.
//   Admitting multi-pattern observables means the tracker now indexes
//   `pat.molecules[0]` for every pattern of every observable it takes,
//   and one pattern shape has no [0]: a `<Pattern>` whose
//   `<ListOfMolecules>` is empty.  `evaluate_observable` has always
//   skipped such a pattern, but the tracker's seeding walk read the
//   empty vector and the process segfaulted before the first sample —
//   already true for a `Molecules` observable before this change, since
//   those were never gated on pattern count.  BNGL cannot express the
//   shape and BNG2 will never write it, so only a host handing
//   RuleMonkey XML directly reaches it; empty_pattern_obs_model.xml is
//   that XML, hand-spliced.  An observable carrying one is now left to
//   the full walk entire, which skips the pattern and counts the rest.

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

void test_multi_pattern_species_matches_molecules_twin(const std::string& species_xml,
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

  check(obs_same, "multi-pattern Species-obs and Molecules-obs arms should give identical Rtot "
                  "trajectories (a rate law reading a stale Species obs separates them)");
  check(fn_same, "multi-pattern Species-obs and Molecules-obs arms should give identical "
                 "conv_rate values");

  if (!obs_same) {
    for (size_t i = 0; i < s_obs.size() && i < m_obs.size(); ++i)
      std::fprintf(stderr, "  row %2zu: species=%.17g molecules=%.17g\n", i, s_obs[i], m_obs[i]);
  }

  // The observable has to actually carry the overlap it is here to test:
  // R() matches every R and R(x~U)/R(x~P) partition the same pool, so a
  // merged-per-complex tally would land on Rtot = R + S rather than
  // 2R + S.  The Molecules twin is the arithmetic being matched.
  check(!s_obs.empty() && s_obs.front() > 0.0,
        "Rtot should be non-zero at t=0, or the trajectory comparison is vacuous");

  // The defect.  Floor is ~1.7; without the fix the Species arm does not
  // finish inside two minutes.
  constexpr double kMaxRatio = 10.0;
  double const ratio = t_species / (t_molecules > 0 ? t_molecules : 1e-9);
  check(ratio < kMaxRatio,
        "a rate law over a multi-pattern Species obs should cost about what the same rate law "
        "over a multi-pattern Molecules obs costs; ratio was " +
            std::to_string(ratio) + " (species " + std::to_string(t_species) + "s, molecules " +
            std::to_string(t_molecules) + "s)");
  std::fprintf(stderr, "arm 1: species %.4fs, molecules %.4fs, ratio %.2f\n", t_species,
               t_molecules, ratio);
}

// ---- ARM 2 -------------------------------------------------------------

void test_multi_pattern_species_cx_bookkeeping(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  double const t_tot = sim.get_parameter("T_tot");
  double const b_tot = sim.get_parameter("B_tot");
  double const kt = sim.get_parameter("kt");
  auto const r = sim.run({0.0, 60.0, 12}, kSeed);

  const auto& mixed = series(r, "Mixed");
  const auto& mixed_walk = series(r, "MixedWalk");
  const auto& chains = series(r, "Chains");
  const auto& chains_ctl = series(r, "ChainsCtl");
  const auto& a_bound = series(r, "A_bound");
  const auto& rate = fn_series(r, "mixed_rate");

  check(!mixed.empty() && mixed.size() == mixed_walk.size() && mixed.size() == chains.size() &&
            mixed.size() == chains_ctl.size(),
        "every observable should have the same number of sample rows");

  // The semantics, pinned on the seed.  Each seeded trimer carries three
  // A's and one B, so it satisfies both of Mixed's patterns and both of
  // Chains'; each free B satisfies Mixed's `B()==1` alone.  A tally
  // merged across an observable's patterns would give t_tot + b_tot and
  // t_tot instead.  This row is the tracker's seeded value, not a walk:
  // refresh_observables_for_sample leaves tracked observables alone.
  double const expect_mixed0 = (2.0 * t_tot) + b_tot;
  double const expect_chains0 = 2.0 * t_tot;
  if (!mixed.empty()) {
    check(mixed.front() == expect_mixed0,
          "Mixed at t=0 should be 2*T_tot + B_tot = " + std::to_string(expect_mixed0) +
              " (one count per pattern the complex satisfies), got " +
              std::to_string(mixed.front()));
    check(chains.front() == expect_chains0,
          "Chains at t=0 should be 2*T_tot = " + std::to_string(expect_chains0) +
              " (one count per pattern the complex satisfies), got " +
              std::to_string(chains.front()));
    check(a_bound.front() == expect_chains0,
          "A_bound at t=0 should be 2*T_tot = " + std::to_string(expect_chains0) + ", got " +
              std::to_string(a_bound.front()));
  }

  // Same instant, from-scratch walk: the oracle the seeded value has to
  // match, taken on a fresh session so nothing has been delta-updated.
  {
    rulemonkey::RuleMonkeySimulator fresh(xml);
    fresh.initialize(kSeed);
    auto const truth = fresh.get_observable_values();
    const auto& names = fresh.observable_names();
    int const i_mixed = idx_of(names, "Mixed");
    int const i_chains = idx_of(names, "Chains");
    check(i_mixed >= 0 && i_chains >= 0, "Mixed and Chains should exist as observables");
    if (i_mixed >= 0)
      check(truth[static_cast<size_t>(i_mixed)] == expect_mixed0,
            "the from-scratch walk of Mixed at t=0 should be 2*T_tot + B_tot too, got " +
                std::to_string(truth[static_cast<size_t>(i_mixed)]));
    if (i_chains >= 0)
      check(truth[static_cast<size_t>(i_chains)] == expect_chains0,
            "the from-scratch walk of Chains at t=0 should be 2*T_tot too, got " +
                std::to_string(truth[static_cast<size_t>(i_chains)]));
  }

  // Each read observable against its unread twin, at every sample.
  // MixedWalk is served by the full walk, so this is tracker vs oracle on
  // every row; ChainsCtl is tracked on the per-sample cadence, so this is
  // the two flush cadences against each other over a shared dead-complex
  // side channel.
  struct Pair {
    const char* read;
    const char* twin;
    const std::vector<double>* a;
    const std::vector<double>* b;
  };
  const Pair pairs[] = {{"Mixed", "MixedWalk", &mixed, &mixed_walk},
                        {"Chains", "ChainsCtl", &chains, &chains_ctl}};
  for (const auto& p : pairs) {
    bool agree = p.a->size() == p.b->size();
    for (size_t i = 0; i < p.a->size() && i < p.b->size(); ++i) {
      if ((*p.a)[i] != (*p.b)[i]) {
        agree = false;
        std::fprintf(stderr, "  row %2zu: %s=%.17g %s=%.17g\n", i, p.read, (*p.a)[i], p.twin,
                     (*p.b)[i]);
      }
    }
    check(agree, std::string("the rate-read '") + p.read + "' (per-event flush) and its unread " +
                     "twin '" + p.twin + "' must hold the same value at every sample");
  }

  // The rate law's own column must be the observables it reads, times kt.
  bool rate_ok = rate.size() == mixed.size();
  for (size_t i = 0; i < rate.size() && i < mixed.size() && i < chains.size(); ++i)
    if (std::fabs(rate[i] - (kt * (mixed[i] + chains[i]))) > 1e-12 * (1.0 + std::fabs(rate[i])))
      rate_ok = false;
  check(rate_ok, "mixed_rate column should equal kt*(Mixed + Chains) at every sample");

  // Independent oracle at t_end: a from-scratch walk of every observable,
  // taken on a live session driven to the same instant with the same
  // seed.  Nothing during a run ever re-derives a tracked value that way,
  // so a delta the tracker got wrong at any event of the run is still
  // there to be caught — and the call is made only after the last step,
  // because it re-bases obs_values and would heal any drift it did not
  // first report.
  //
  // Every observable is checked, not just the Species ones: A_bound and
  // Q_tot ride the long-settled Molecules delta path, so if they agree
  // and the Species ones do not, the tracker is at fault rather than the
  // trajectory.
  rulemonkey::RuleMonkeySimulator live(xml);
  live.initialize(kSeed);
  live.step_to(60.0);
  auto const truth = live.get_observable_values();
  const auto& names = live.observable_names();
  for (const char* nm : {"Mixed", "MixedWalk", "Chains", "ChainsCtl", "A_bound", "Q_tot"}) {
    int const i = idx_of(names, nm);
    check(i >= 0, std::string("observable '") + nm + "' should exist");
    if (i < 0 || mixed.empty())
      continue;
    double const tracked = series(r, nm).back();
    check(truth[static_cast<size_t>(i)] == tracked,
          std::string("incrementally tracked '") + nm + "' (" + std::to_string(tracked) +
              ") should equal a from-scratch walk (" +
              std::to_string(truth[static_cast<size_t>(i)]) + ") at t_end");
  }

  // Sanity: the model has to actually exercise complexes and the rule
  // the observables price, or the assertions above are vacuous.  The
  // seed already forms chains, so the check that matters is that the
  // chemistry moved them rather than that they exist.
  check(!chains.empty() && chains.back() > 0, "the model should still hold chains at t_end");
  check(!a_bound.empty() && a_bound.back() > a_bound.front(),
        "polymerisation should have bound more A than the seed did");
  check(series(r, "Q_tot").back() > 0, "the rate-priced rule should have fired");
}

// ---- ARM 3 -------------------------------------------------------------

void test_empty_pattern_observable(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  auto const r = sim.run({0.0, 5.0, 5}, kSeed);

  const auto& rmix = series(r, "Rmix");   // empty pattern + R()
  const auto& rspec = series(r, "Rspec"); // empty pattern + R(), Species
  const auto& rnone = series(r, "Rnone"); // the empty pattern alone
  const auto& rtot = series(r, "Rtot");   // R() alone

  check(!rtot.empty() && rtot.front() > 0, "the model should start with a non-empty R pool");

  // An empty pattern contributes nothing, so pairing one with `R()`
  // must leave the observable equal to `R()` alone, and an observable
  // that is nothing but an empty pattern must stay at zero.
  bool mix_ok = rmix.size() == rtot.size();
  bool spec_ok = rspec.size() == rtot.size();
  bool none_ok = rnone.size() == rtot.size();
  for (size_t i = 0; i < rtot.size(); ++i) {
    if (i >= rmix.size() || rmix[i] != rtot[i])
      mix_ok = false;
    if (i >= rspec.size() || rspec[i] != rtot[i])
      spec_ok = false;
    if (i >= rnone.size() || rnone[i] != 0.0)
      none_ok = false;
  }
  check(mix_ok, "a Molecules observable pairing an empty pattern with R() should equal R() alone");
  check(spec_ok, "a Species observable pairing an empty pattern with R() should equal R() alone "
                 "(every R is a free monomer)");
  check(none_ok, "an observable whose only pattern is empty should stay at zero");

  // The rate law reads Rmix, so its column is the proof that the value
  // the propensity saw is the same one the sample row reports.
  double const kt = sim.get_parameter("kt");
  const auto& rate = fn_series(r, "conv_rate");
  bool rate_ok = rate.size() == rmix.size();
  for (size_t i = 0; i < rate.size() && i < rmix.size(); ++i)
    if (std::fabs(rate[i] - (kt * rmix[i])) > 1e-12 * (1.0 + std::fabs(rate[i])))
      rate_ok = false;
  check(rate_ok, "conv_rate column should equal kt*Rmix at every sample");

  // And the from-scratch walk, which is what the fallback path must be
  // agreeing with.
  rulemonkey::RuleMonkeySimulator live(xml);
  live.initialize(kSeed);
  live.step_to(5.0);
  auto const truth = live.get_observable_values();
  const auto& names = live.observable_names();
  for (const char* nm : {"Rmix", "Rspec", "Rnone", "Rtot"}) {
    int const i = idx_of(names, nm);
    check(i >= 0, std::string("observable '") + nm + "' should exist");
    if (i < 0 || rtot.empty())
      continue;
    check(truth[static_cast<size_t>(i)] == series(r, nm).back(),
          std::string("'") + nm + "' (" + std::to_string(series(r, nm).back()) +
              ") should equal a from-scratch walk (" +
              std::to_string(truth[static_cast<size_t>(i)]) + ") at t_end");
  }
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 5) {
    std::fprintf(
        stderr,
        "usage: %s <species_obs.xml> <molecules_obs.xml> <species_cx.xml> <empty_pattern.xml>\n",
        argv[0]);
    return 2;
  }
  try {
    test_multi_pattern_species_matches_molecules_twin(argv[1], argv[2]);
    test_multi_pattern_species_cx_bookkeeping(argv[3]);
    test_empty_pattern_observable(argv[4]);
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
