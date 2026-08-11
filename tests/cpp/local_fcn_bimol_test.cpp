// DOR1 rate-parity test: a bimolecular rule whose rate is a single local
// function, with only one of the two reactants tagged (issue #34).
//
// Before the fix, `recompute_rule_state` had branches for FunctionProduct
// (two tagged reactants) and for `is_local` with molecularity <= 1, but none
// for `is_local` with molecularity 2.  Such a rule fell through to the
// mass-action path, which left `rs.local_propensity_total` at zero while
// `rs.has_local_rates` stayed true — so the first `incremental_update` read
// that never-populated accumulator as the rule's propensity.  The rule fired
// once (off the mass-action propensity computed at load) and then went inert,
// with no clamp warning and a conserved, entirely plausible-looking
// trajectory.
//
// Model: tests/models/feature_coverage/nf_local_fcn_bimol.bngl, four
// independent substrates driven by the four shapes a single tag can take:
//
//   RB  Sa(s~0) + E()%x                 tag on pattern 1, molecule scope
//   RA  E()%x + Sb(s~0)                 tag on pattern 0, molecule scope
//   RC  Sc(s~0) + %x:E()                tag on pattern 1, complex scope
//   RD  Sd(s~0) + %x:E(d!1).E(d!1)      complex scope, multi-molecule pattern
//
// The enzyme pool is inert (no rule touches E), so each substrate sees a
// constant per-molecule hazard c = Σ_t f(t) over the tagged reactant's
// matches, and Sub(t) is an exact binomial death process:
//
//     Sub(t) ~ Binomial(S0, exp(-c·t))
//
// With kc = 1e-3, 20 free E(m~0) and 10 E(m~1) homodimers (Emod = 20):
//
//   molecule scope  c = kc · Σ_{40 E molecules} Emod(x)  = kc·20 = 0.02
//   complex scope   c = kc · Σ_{20 dimer E's}   Emod(cx) = kc·40 = 0.04
//
// and for RD the same 0.04 arrives as 10 dimers × 2 embeddings × 2.  That
// gives an analytic mean and variance to test against, with no NFsim or BNG2
// dependency — the harness model nf_local_fcn_bimol covers the NFsim
// ensemble comparison separately.

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
  int const i = idx_of(r, name);
  if (i < 0)
    throw std::runtime_error("observable not found: " + name);
  return r.observable_data[i].back();
}

// One tested rule: its substrate observable and the per-molecule hazard the
// DOR1 propensity is supposed to realize.
struct Arm {
  const char* observable;
  double hazard;
  const char* shape;
};

constexpr double kS0 = 300.0;  // seed count of each substrate
constexpr double kTEnd = 60.0; // matches the model's simulate action
constexpr int kReps = 60;      // replicates per arm
constexpr double kSigmaTol = 4.0;

void test_dor1_arms(const std::string& xml) {
  const std::vector<Arm> arms = {
      {"SubA", 0.02, "tag on reactant 1, molecule scope"},
      {"SubB", 0.02, "tag on reactant 0, molecule scope"},
      {"SubC", 0.04, "tag on reactant 1, complex scope"},
      {"SubD", 0.04, "complex scope, multi-molecule tagged pattern"},
  };

  std::vector<double> sum(arms.size(), 0.0);
  int stalled = 0;

  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    auto r = sim.run({0.0, kTEnd, 2},
                     /*seed=*/std::uint64_t{7000} + static_cast<std::uint64_t>(rep));

    // The enzyme pool is context only — no rule may consume or modify it.
    check(final_value(r, "Etot") == 40.0, "Etot conserved at 40");
    check(final_value(r, "Emod") == 20.0, "Emod conserved at 20");

    for (size_t a = 0; a < arms.size(); ++a) {
      double const v = final_value(r, arms[a].observable);
      sum[a] += v;
      // The pre-fix signature, named explicitly: a single firing and then a
      // dead rule.  Every arm's expected terminal value is <= 91 of 300, so
      // anything above 295 is the rule having stopped, not a tail draw.
      if (v > 295.0)
        ++stalled;
    }
  }

  check(stalled == 0, "no arm stalled after its first firing");

  for (size_t a = 0; a < arms.size(); ++a) {
    double const p = std::exp(-arms[a].hazard * kTEnd);
    double const analytic = kS0 * p;
    double const se = std::sqrt(kS0 * p * (1.0 - p) / kReps);
    double const empirical = sum[a] / kReps;
    double const diff = std::fabs(empirical - analytic);

    std::fprintf(stderr, "%-5s (%-44s): empirical=%7.3f analytic=%7.3f |diff|=%6.3f (%.2f SE)\n",
                 arms[a].observable, arms[a].shape, empirical, analytic, diff, diff / se);
    check(diff < kSigmaTol * se,
          std::string(arms[a].observable) + " terminal mean within 4 SE of the binomial analytic");
  }
}

// The two scopes must not collapse onto each other: a complex-wide arm has to
// run at exactly twice the molecule-scope hazard, which is what distinguishes
// "the tag was resolved to the right scope" from "the tag was ignored".
void test_scope_separation(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  auto r = sim.run({0.0, kTEnd, 2}, /*seed=*/99);
  double const mol_scope = final_value(r, "SubA");
  double const cx_scope = final_value(r, "SubC");
  std::fprintf(stderr, "scope separation: SubA=%.0f SubC=%.0f\n", mol_scope, cx_scope);
  check(cx_scope < mol_scope, "complex-wide arm decays faster than the molecule-scope arm");
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "Usage: local_fcn_bimol_test <nf_local_fcn_bimol.xml>\n");
    return 2;
  }
  std::string const xml = argv[1];
  try {
    test_dor1_arms(xml);
    test_scope_separation(xml);
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
