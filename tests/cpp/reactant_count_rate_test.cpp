// `reactant_N()` in a rate law (issue #59).
//
// A BNGL model that wants the match count of a rule's Nth reactant pattern
// inside that rule's own rate law declares an empty function by that name
// and writes it into the rate expression:
//
//     begin functions
//        reactant_1()
//        reactant_2()
//        NucF()=if(Dimer<1,kNuc,0)
//     end functions
//     A(b) + A(b) -> A(b!1).A(b!1)   reactant_1()*reactant_2()*NucF()
//
// BNG2 emits both placeholders as ordinary functions with an empty
// <Expression>, so the XML carries no hint of what they mean; the name is
// the whole convention, and NFsim reads it the same way.  RM used to read
// them as what they literally are — a function with no body, i.e. the
// constant zero — so every rate law built on one evaluated to zero.  The
// rule never fired, the run ended with no events, and nothing was reported:
// a fast, clean-looking, entirely wrong trajectory.
//
// The rate law is not a total rate, so the count the placeholder supplies is
// multiplied in ON TOP of the ordinary mass-action count.  For the
// nucleation model above, with 200 free A and kNuc = 0.01:
//
//     rate      = reactant_1 * reactant_2 * NucF = 200 * 200 * 0.01 = 400
//     propensity = (200*200 - 200)/2 * 400       = 7.96e6
//
// (the homodimer count is the number of distinct unordered pairs, which is
// what compute_propensity's same_components branch realizes).
//
// Three things are checked here:
//
//   1. The reported value of the rate law is that 400 exactly, and the
//      placeholders themselves report zero — they have no value outside the
//      rule that asks for them.
//   2. The PROPENSITY is 7.96e6 and not merely non-zero.  NucF gates on
//      Dimer, so the rule fires at most once whatever its rate; running many
//      replicates to exactly the mean waiting time 1/a and counting how many
//      of them fired turns that rate into an observable number, which must
//      come out at 1 - 1/e.
//   3. On the feature-coverage model, each arm written with placeholders
//      tracks a control arm that states the same law through ordinary
//      observables, and both match the closed-form solution.
//
// Plus the two shapes RM refuses at load rather than resolving to zero: a
// count of a reactant the rule does not have, and a count read from a rate
// law that is also a local (per-instance) function.

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

double final_observable(const rulemonkey::Result& r, const std::string& name) {
  for (size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == name)
      return r.observable_data[i].back();
  throw std::runtime_error("observable not found: " + name);
}

double first_function_value(const rulemonkey::Result& r, const std::string& name) {
  for (size_t i = 0; i < r.function_names.size(); ++i)
    if (r.function_names[i] == name)
      return r.function_data[i].front();
  throw std::runtime_error("function not found: " + name);
}

// ---------------------------------------------------------------------------
// 1 + 2: the reported rate law, and the propensity behind it
// ---------------------------------------------------------------------------

constexpr double kNuc = 0.01;
constexpr double kA0 = 200.0;
// reactant_1 * reactant_2 * NucF at the seed state.
constexpr double kExpectedRate = kA0 * kA0 * kNuc;
// Distinct unordered pairs of free A, times that rate.
constexpr double kExpectedPropensity = ((kA0 * kA0) - kA0) / 2.0 * kExpectedRate;

void test_nucleation_rate_value(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  auto r = sim.run({0.0, 5.0, 11}, /*seed=*/1);

  double const rate = first_function_value(r, "_rateLaw1");
  std::fprintf(stderr, "_rateLaw1 at t=0: %.6f (expected %.6f)\n", rate, kExpectedRate);
  check(std::fabs(rate - kExpectedRate) < 1e-9, "_rateLaw1 reports reactant_1*reactant_2*NucF");

  // The placeholders have no value of their own: they stand for whichever
  // rule is asking, so the settled function table leaves them at zero.
  check(first_function_value(r, "reactant_1") == 0.0, "reactant_1 reports no global value");
  check(first_function_value(r, "reactant_2") == 0.0, "reactant_2 reports no global value");

  // The pre-fix signature: nothing ever fires.  NucF closes the gate after
  // the first dimer, so Dimer counts 2 molecules and stops there.
  check(final_observable(r, "Dimer") == 2.0, "one dimer forms and the gate closes");
  check(final_observable(r, "Atot") == kA0, "no A is created or destroyed");
}

// Run each replicate for exactly the mean waiting time of the expected
// propensity.  The rule can fire at most once, so the fraction of
// replicates that fired estimates 1 - exp(-a*t) = 1 - 1/e, and the estimate
// depends on the propensity and on nothing else.
void test_nucleation_propensity(const std::string& xml) {
  constexpr int kReps = 400;
  double const t_end = 1.0 / kExpectedPropensity;
  double const p_expect = 1.0 - std::exp(-1.0);
  double const se = std::sqrt(p_expect * (1.0 - p_expect) / kReps);

  int fired = 0;
  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    auto r = sim.run({0.0, t_end, 2},
                     /*seed=*/std::uint64_t{5900} + static_cast<std::uint64_t>(rep));
    if (final_observable(r, "Dimer") > 0.0)
      ++fired;
  }

  double const p = static_cast<double>(fired) / kReps;
  std::fprintf(stderr, "fired in %d/%d replicates at t = 1/a = %.3e: p=%.4f (expected %.4f, %.2f SE)\n",
               fired, kReps, t_end, p, p_expect, std::fabs(p - p_expect) / se);
  check(std::fabs(p - p_expect) < 4.0 * se,
        "firing fraction at the mean waiting time matches the expected propensity");
}

// ---------------------------------------------------------------------------
// 3: placeholder arms against control arms on the coverage model
// ---------------------------------------------------------------------------

constexpr double kCovTEnd = 100.0;
constexpr int kCovReps = 40;

void test_coverage_arms(const std::string& xml) {
  // Closed forms for the two laws, both of which the placeholder arm and the
  // control arm state (see ft_reactant_count_rate.bngl).
  constexpr double ku = 2.0e-4;
  constexpr double kb = 2.5e-8;
  double const uni_expect = 200.0 / (1.0 + (ku * 200.0 * kCovTEnd));
  double const bim_expect = std::pow(std::pow(100.0, -3.0) + (3.0 * kb * kCovTEnd), -1.0 / 3.0);

  double u_sum = 0, v_sum = 0, a_sum = 0, c_sum = 0;
  for (int rep = 0; rep < kCovReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    auto r = sim.run({0.0, kCovTEnd, 2},
                     /*seed=*/std::uint64_t{5910} + static_cast<std::uint64_t>(rep));
    u_sum += final_observable(r, "Uoff");
    v_sum += final_observable(r, "Voff");
    a_sum += final_observable(r, "Afree");
    c_sum += final_observable(r, "Cfree");
  }

  double const u = u_sum / kCovReps;
  double const v = v_sum / kCovReps;
  double const a = a_sum / kCovReps;
  double const c = c_sum / kCovReps;
  std::fprintf(stderr, "unimolecular: reactant_1() arm=%.2f  observable arm=%.2f  closed form=%.2f\n",
               u, v, uni_expect);
  std::fprintf(stderr, "bimolecular:  reactant_N() arm=%.2f  observable arm=%.2f  closed form=%.2f\n",
               a, c, bim_expect);

  // A count applied zero times leaves the arm at its seed value; a count
  // applied once instead of twice leaves it near it (196 of 200, 100 of 100).
  // 5% of the closed form is far tighter than either.
  check(std::fabs(u - uni_expect) < 0.05 * uni_expect,
        "unimolecular reactant_1() arm matches the closed form");
  check(std::fabs(a - bim_expect) < 0.05 * bim_expect,
        "bimolecular reactant_1()*reactant_2() arm matches the closed form");
  check(std::fabs(u - v) < 0.05 * uni_expect, "unimolecular arm tracks its observable control");
  check(std::fabs(a - c) < 0.05 * bim_expect, "bimolecular arm tracks its observable control");
}

// ---------------------------------------------------------------------------
// The two shapes RM refuses rather than resolving to zero
// ---------------------------------------------------------------------------

void test_refusal(const std::string& xml, const std::string& must_contain,
                  const std::string& tag) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  bool found = false;
  for (const auto& u : sim.unsupported_features()) {
    if (u.severity == rulemonkey::Severity::Error &&
        u.feature.find(must_contain) != std::string::npos)
      found = true;
  }
  if (!found) {
    std::fprintf(stderr, "  %s: reported features were:\n", tag.c_str());
    for (const auto& u : sim.unsupported_features())
      std::fprintf(stderr, "    [%d] %s\n", static_cast<int>(u.severity), u.feature.c_str());
  }
  check(found, tag + " is refused at load with an Error");
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 5) {
    std::fprintf(stderr, "Usage: reactant_count_rate_test <nucleation.xml> <coverage.xml> "
                         "<arity.xml> <local_fcn.xml>\n");
    return 2;
  }
  try {
    test_nucleation_rate_value(argv[1]);
    test_nucleation_propensity(argv[1]);
    test_coverage_arms(argv[2]);
    test_refusal(argv[3], "reactant_2()", "a count of a reactant the rule does not have");
    test_refusal(argv[4], "local (per-instance) function",
                 "a count read from a local-function rate law");
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
