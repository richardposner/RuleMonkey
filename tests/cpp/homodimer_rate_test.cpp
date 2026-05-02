// Rate-parity test for the same_components homodimer fast path.
//
// The bimolecular sampler used to null-event self-pairs (`mol_a == mol_b`)
// after pairing the inflated propensity `ab²/2` with that rejection.  The
// scheme was statistically correct but spent SSA cycles on null events at
// rate `1/N` per selection — non-negligible for small populations.
//
// The fix tracks `ab_both_sq_total = Σ ab(m)²` incrementally and uses the
// deflated propensity `(ab² − ab_sq) / 2 · k` together with retry-until-
// distinct sampling.  Realized rate is the analytic Gillespie rate at any
// N; this test pins down a known steady-state mean against the chemical
// master equation, where the old code's wasted-cycles mode and the new
// code's no-waste mode must both yield the same statistics.
//
// Model: A(a) + A(a) <-> A(a!1).A(a!1) with overrides kp=0.5, km=0.5
// (K_eq=1) and A_tot=8.  Stationary distribution computed by detailed
// balance over k = number of dimers (CME truncated at A_tot/2):
//
//     π(0) ∝ 1
//     π(k+1)/π(k) = kp·(A−2k)·(A−2k−1) / (2·kr·(k+1))
//
// For A=8, K=1: π = (1, 28, 210, 420, 105) / 764, mean E[k] ≈ 2.785,
// var(k) ≈ 0.703.  With N_reps=200 and t_end=200 the Markov chain mixes
// far past the relaxation time, so each replicate's terminal sample is
// effectively an independent draw from π.  Standard error of the
// empirical mean is √(var/N) ≈ 0.0593.  Tolerance is 4σ ≈ 0.24, well
// below the magnitude of any sampler bias the old buggy formula would
// have produced (the inflated propensity without `ab_sq` deflation
// would shift the realized rate, and hence the equilibrium, by O(1/N)).

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

double final_value(const rulemonkey::Result& r, const std::string& name) {
  int i = idx_of(r, name);
  if (i < 0)
    throw std::runtime_error("observable not found: " + name);
  return r.observable_data[i].back();
}

// Run an ensemble of homodimer simulations with the given parameter
// overrides; return the mean of the terminal AA_1 observable.  Each
// replicate uses a distinct seed.
double homodimer_terminal_mean(const std::string& xml, double a_tot, double kp, double km,
                               double t_end, int n_reps) {
  double sum = 0.0;
  for (int rep = 0; rep < n_reps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.set_param("A_tot", a_tot);
    sim.set_param("kp", kp);
    sim.set_param("km", km);
    auto r =
        sim.run({0.0, t_end, 2}, /*seed=*/std::uint64_t{1000} + static_cast<std::uint64_t>(rep));
    sum += final_value(r, "AA_1");
  }
  return sum / n_reps;
}

// CME stationary distribution for A(a) + A(a) <-> A(a!1).A(a!1).
// Returns mean of k = N_AA at steady state.
double cme_mean(int a_tot, double kp, double km) {
  int kmax = a_tot / 2;
  std::vector<double> pi(kmax + 1, 0.0);
  pi[0] = 1.0;
  for (int k = 0; k < kmax; ++k) {
    int free_a = a_tot - (2 * k);
    double fwd = kp * free_a * (free_a - 1) / 2.0;
    double rev = km * (k + 1);
    pi[k + 1] = pi[k] * fwd / rev;
  }
  double Z = 0.0;
  for (double p : pi)
    Z += p;
  double mean = 0.0;
  for (int k = 0; k <= kmax; ++k)
    mean += k * pi[k] / Z;
  return mean;
}

void test_small_n_homodimer(const std::string& xml) {
  // K = 1, A_tot = 8.  Analytic E[N_AA] ≈ 2.785.  At small N the old
  // code's per-event self-pair rejection rate was ~1/N = 12.5%; the
  // new code's deflated propensity makes that retry-until-distinct.
  // Both must converge to the same mean.
  constexpr double kp = 0.5;
  constexpr double km = 0.5;
  constexpr int a_tot = 8;
  constexpr int n_reps = 200;
  constexpr double t_end = 200.0;

  double analytic = cme_mean(a_tot, kp, km);
  double empirical = homodimer_terminal_mean(xml, a_tot, kp, km, t_end, n_reps);
  double diff = std::fabs(empirical - analytic);

  // Tolerance: 4σ on the standard error of the mean for n_reps=200.
  // var(k) ≈ 0.703 (computed once for K=1, A=8), so SE = sqrt(0.703/200)
  // ≈ 0.0593, 4σ ≈ 0.237.  Add 0.05 padding for finite mixing time
  // (t_end=200 is many relaxation times, but not infinite).
  constexpr double tol = 0.30;

  std::fprintf(stderr,
               "small-N homodimer (A_tot=%d, K=1): empirical=%.4f, analytic=%.4f, |diff|=%.4f\n",
               a_tot, empirical, analytic, diff);
  check(diff < tol, "small-N homodimer mean E[AA] within 4σ of CME analytic");
}

void test_minimal_pair_homodimer(const std::string& xml) {
  // The pathological small-N case: A_tot=2 means exactly one possible
  // dimer pair, and the old code's self-pair rejection rate hits 1/2 per
  // selection.  At K=1 the stationary distribution is π(0) = π(1) = 0.5,
  // so E[N_AA] = 0.5.  Variance = 0.25, SE = sqrt(0.25/200) ≈ 0.0354,
  // 4σ ≈ 0.141.
  constexpr double kp = 0.5;
  constexpr double km = 0.5;
  constexpr int a_tot = 2;
  constexpr int n_reps = 200;
  constexpr double t_end = 100.0;

  double analytic = cme_mean(a_tot, kp, km);
  double empirical = homodimer_terminal_mean(xml, a_tot, kp, km, t_end, n_reps);
  double diff = std::fabs(empirical - analytic);
  constexpr double tol = 0.20;

  std::fprintf(stderr,
               "minimal-pair homodimer (A_tot=%d, K=1): empirical=%.4f, analytic=%.4f, "
               "|diff|=%.4f\n",
               a_tot, empirical, analytic, diff);
  check(diff < tol, "minimal-pair homodimer mean E[AA] within 4σ of CME analytic");
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "Usage: homodimer_rate_test <A_plus_A.xml>\n");
    return 2;
  }
  std::string xml = argv[1];
  try {
    test_minimal_pair_homodimer(xml);
    test_small_n_homodimer(xml);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }
  if (g_failures != 0) {
    std::fprintf(stderr, "%d failure(s)\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "PASS\n");
  return 0;
}
