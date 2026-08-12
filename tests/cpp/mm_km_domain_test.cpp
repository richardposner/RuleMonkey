// The Km <= 0 boundary of the MM(kcat,Km) rate law (issue #46).
//
// The branch guarded its division with
//
//     if (Km + sFree <= 0) return 0;
//
// which answers 0 where the limit is not 0, and does not answer at all where
// the expression is NaN.
//
// Km == 0 is a removable singularity.  sFree = max(S-E, 0) there, so S > E
// already reads kcat*E, but S <= E reads 0/0 — and
//
//     sFree -> Km*S/(E-S)  as Km -> 0+   =>   a -> kcat*S
//
// so the law on the whole Km = 0 line is a = kcat*min(S,E): binding
// infinitely tight, turnover substrate-limited when the enzyme is in excess
// and enzyme-limited otherwise.  Returning 0 made the rule silently inert,
// and froze a model that started at S > E as soon as its substrate decayed
// to S == E — the trajectory simply stopped, mass conserved, no warning.
//
// The approach to that limit was also wrong before it was reached, because
// 0.5*(diff + sqrt(diff^2 + 4*Km*S)) cancels catastrophically for diff < 0.
// On S=100, E=200, kcat=1, where the limit is 100, the old form read
// 100.093 at Km=1e-12, 117.392 at Km=1e-14, and 0 from about 1e-15 down.
// The conjugate form 2*Km*S/(q - diff) is algebraically identical and holds
// the limit exactly.  Km is a parameter, so a scan or a fit walks through
// that band on its way to zero.
//
// Km < 0 is outside the law entirely.  The old guard let the resulting NaN
// through as the rule's propensity — `NaN < 0` is false and std::max
// returned it unchanged — so total_propensity became NaN and the sampler
// kept drawing against it: an MM rule with Km = -1 fired 80 events, exited
// 0, and wrote an ordinary-looking .gdat.  A statically negative Km is now
// refused at load; one arriving later through an override is clamped to zero
// with one warning naming Km.
//
// Model: tests/models/feature_coverage/xml/ft_mm_ratelaw.xml, an ordinary
// asymmetric MM rule (S + E -> P + E), driven entirely through set_param.
// Oracles are closed form, no external simulator:
//
//   Km = 0, S < E   a = kcat*S, so each substrate molecule dies independently
//                   at hazard kcat and S(t) ~ Binomial(S0, exp(-kcat*t)) —
//                   mean and variance both exact.
//   Km = 0, S > E   a = kcat*E while S > E (mean exactly S0 - kcat*E*t) and
//                   kcat*S after, both linear in the count, so the mean-field
//                   integration below is exact away from the crossover.
//   Km = 1e-14/1e-16, S < E  same exponential: the limit is approached, not
//                   jumped to.  The old form reads 17% fast at 1e-14 and
//                   dead at 1e-16.

#include "rulemonkey/simulator.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
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

int idx_of(const rulemonkey::Result& r, const std::string& name) {
  for (size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == name)
      return static_cast<int>(i);
  return -1;
}

double value_at(const rulemonkey::Result& r, const std::string& name, size_t t_idx) {
  int const i = idx_of(r, name);
  if (i < 0)
    throw std::runtime_error("observable not found: " + name);
  return r.observable_data[i][t_idx];
}

// Mean-field substrate trajectory under a = kcat*min(S,E), the Km = 0 law.
double expected_substrate(double s0, double e0, double kcat, double T) {
  constexpr int kSteps = 200000;
  double const dt = T / kSteps;
  double s = s0;
  for (int i = 0; i < kSteps; ++i)
    s -= kcat * std::min(s, e0) * dt;
  return s;
}

struct Ensemble {
  double mean = 0;
  double se = 0;
};

Ensemble run_ensemble(rulemonkey::RuleMonkeySimulator& sim, const rulemonkey::TimeSpec& ts,
                      const std::string& obs, size_t t_idx, int reps) {
  double sum = 0;
  double sumsq = 0;
  for (int rep = 0; rep < reps; ++rep) {
    rulemonkey::Result const r = sim.run(ts, std::uint64_t{4000} + static_cast<std::uint64_t>(rep));
    double const v = value_at(r, obs, t_idx);
    sum += v;
    sumsq += v * v;
  }
  double const mean = sum / reps;
  double const var = (sumsq / reps) - (mean * mean);
  return {mean, std::sqrt((var > 0 ? var : 0.0) / reps)};
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <ft_mm_ratelaw.xml>\n", argv[0]);
    return 2;
  }
  std::string const xml = argv[1];
  constexpr int kReps = 200;
  constexpr double kKcat = 1.0;

  rulemonkey::TimeSpec ts;
  ts.t_start = 0.0;
  ts.t_end = 3.0;
  ts.n_points = 3; // t = 0, 1, 2, 3

  try {
    rulemonkey::RuleMonkeySimulator sim(xml);

    // --- Km = 0 with the enzyme in excess: a = kcat*S, a pure death process.
    // The old code returned 0 here and the rule never fired at all.
    for (double const km : {0.0, 1e-14, 1e-16}) {
      sim.clear_param_overrides();
      sim.set_param("S0", 100);
      sim.set_param("E0", 200);
      sim.set_param("kcat", kKcat);
      sim.set_param("Km", km);

      // Guard the driver itself: the overrides must have reached the seeds.
      rulemonkey::Result const probe = sim.run(ts, 1);
      check(value_at(probe, "S_n", 0) == 100.0,
            "set_param(S0) did not reach the seed species (t=0 S_n = " +
                std::to_string(value_at(probe, "S_n", 0)) + ")");

      for (size_t ti = 1; ti <= 3; ++ti) {
        auto const t = static_cast<double>(ti);
        Ensemble const e = run_ensemble(sim, ts, "S_n", ti, kReps);
        // S(t) ~ Binomial(100, exp(-kcat*t)) exactly.
        double const p = std::exp(-kKcat * t);
        double const want = 100.0 * p;
        double const sd = std::sqrt(100.0 * p * (1.0 - p));
        double const tol = (4.0 * sd / std::sqrt(static_cast<double>(kReps))) + (0.01 * want);
        std::printf("Km=%-6g S<E  t=%.0f  S_n=%7.3f  expected=%7.3f  tol=%5.3f\n", km, t, e.mean,
                    want, tol);
        check(std::fabs(e.mean - want) <= tol,
              "Km=" + std::to_string(km) + " S<E at t=" + std::to_string(t) + ": S_n " +
                  std::to_string(e.mean) + " outside " + std::to_string(want) + " +/- " +
                  std::to_string(tol));
      }
    }

    // --- Km = 0 with the substrate in excess: a = kcat*E until S == E, then
    // kcat*S.  The old code ran the first phase correctly and then froze the
    // trajectory at S == E forever.
    {
      sim.clear_param_overrides();
      sim.set_param("S0", 200);
      sim.set_param("E0", 100);
      sim.set_param("kcat", kKcat);
      sim.set_param("Km", 0.0);
      for (size_t ti = 1; ti <= 3; ++ti) {
        auto const t = static_cast<double>(ti);
        Ensemble const e = run_ensemble(sim, ts, "S_n", ti, kReps);
        double const want = expected_substrate(200.0, 100.0, kKcat, t);
        double const tol = (4.0 * e.se) + (0.03 * want);
        std::printf("Km=0      S>E  t=%.0f  S_n=%7.3f  expected=%7.3f  tol=%5.3f\n", t, e.mean,
                    want, tol);
        check(std::fabs(e.mean - want) <= tol,
              "Km=0 S>E at t=" + std::to_string(t) + ": S_n " + std::to_string(e.mean) +
                  " outside " + std::to_string(want) + " +/- " + std::to_string(tol));
      }
      // The freeze regression, stated directly: the substrate must go well
      // below the enzyme count rather than stopping at it.
      Ensemble const tail = run_ensemble(sim, ts, "S_n", 3, kReps);
      check(tail.mean < 30.0, "Km=0 S>E: substrate stalled at S_n=" + std::to_string(tail.mean) +
                                  " (the pre-fix code froze it at S == E = 100)");
    }

    // --- Km < 0 arriving after load: clamped to zero, rule inert, no NaN.
    {
      sim.clear_param_overrides();
      sim.set_param("S0", 200);
      sim.set_param("E0", 100);
      sim.set_param("kcat", kKcat);
      sim.set_param("Km", -1.0);
      rulemonkey::Result const r = sim.run(ts, 7);
      for (size_t ti = 0; ti <= 3; ++ti) {
        double const s = value_at(r, "S_n", ti);
        double const p = value_at(r, "P_n", ti);
        check(!std::isnan(s) && !std::isnan(p),
              "Km<0: NaN reached the observables at t index " + std::to_string(ti));
        check(s == 200.0 && p == 0.0, "Km<0: rule fired (S_n=" + std::to_string(s) +
                                          ", P_n=" + std::to_string(p) +
                                          ") — a negative Km must clamp the propensity to zero");
      }
      std::printf("Km=-1          S_n=%.3f P_n=%.3f at t_end (want 200/0, no NaN)\n",
                  value_at(r, "S_n", 3), value_at(r, "P_n", 3));
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: threw: %s\n", e.what());
    return 1;
  }

  if (g_failures) {
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
  }
  std::printf("mm_km_domain_test: all checks passed\n");
  return 0;
}
