// A rule whose right-hand side can never be true is held inert and named in a
// warning, not refused (issue #105, and its mirror image).
//
// The fixture is the BNG2-oracle model tests/bng_oracle/xml/
// products_never_form.xml, whose counts are pinned against BNG2 there.  This
// test pins what that comparison cannot see:
//
//   * the load reports each such rule as a Warn, so rm_driver runs the model
//     as BNG2 does instead of refusing it, and reports no other rule;
//   * `Join`, `Bare` and `Lone`, which write two molecules nothing joins
//     into one product pattern, and `Apart`, which writes two molecules of
//     one complex into different product patterns and never parts them,
//     never fire under either setting of block_same_complex_binding.  That
//     flag decides only whether one complex may supply both reactant
//     patterns; the `+` between product patterns holds either way (issue
//     #123).

#include "rulemonkey/simulator.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace {

int g_failures = 0;

void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

double final_value(const rulemonkey::Result& r, const std::string& name) {
  for (size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == name)
      return r.observable_data[i].back();
  throw std::runtime_error("observable not found: " + name);
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "Usage: products_never_form_test <products_never_form.xml>\n");
    return 2;
  }

  // k t = 0.5 for the rules that fire at k, and k2 N t = 0.5 for Bound.
  constexpr double kN = 1000.0;
  constexpr double kTEnd = 500.0;
  double const p = 1.0 - std::exp(-0.5);
  double const expected = kN * p;
  double const tol = 5.0 * std::sqrt(kN * p * (1.0 - p));

  try {
    {
      rulemonkey::RuleMonkeySimulator sim(argv[1]);
      int n_warn = 0;
      bool join = false;
      bool bare = false;
      bool lone = false;
      bool apart = false;
      for (const auto& w : sim.unsupported_features()) {
        check(w.severity == rulemonkey::Severity::Warn,
              "a rule that never fires must be a Warn, not an Error: " + w.feature);
        check(w.element == "ReactionRule", "unexpected element " + w.element);
        ++n_warn;
        join = join || w.feature.find("(Join)") != std::string::npos;
        bare = bare || w.feature.find("(Bare)") != std::string::npos;
        lone = lone || w.feature.find("(Lone)") != std::string::npos;
        apart = apart || w.feature.find("(Apart)") != std::string::npos;
      }
      check(n_warn == 4, "expected 4 warnings, got " + std::to_string(n_warn));
      check(join && bare && lone && apart, "the warnings must name Join, Bare, Lone and Apart");

      auto r = sim.run(rulemonkey::TimeSpec{0.0, kTEnd, 2}, 7);
      check(final_value(r, "Join_on") == 0, "Join fired");
      check(final_value(r, "Bare_on") == 0, "Bare fired");
      check(final_value(r, "Lone_on") == 0, "Lone fired");
      check(final_value(r, "Apart_on") == 0, "Apart fired under block_same_complex_binding");
      double const bound = final_value(r, "Bound_on");
      check(std::abs(bound - expected) < tol,
            "Bound (control) should fire: " + std::to_string(bound));
    }
    {
      rulemonkey::RuleMonkeySimulator sim(argv[1]);
      sim.set_block_same_complex_binding(false);
      auto r = sim.run(rulemonkey::TimeSpec{0.0, kTEnd, 2}, 7);
      check(final_value(r, "Join_on") == 0, "Join fired without block_same_complex_binding");
      check(final_value(r, "Bare_on") == 0, "Bare fired without block_same_complex_binding");
      check(final_value(r, "Lone_on") == 0, "Lone fired without block_same_complex_binding");
      check(final_value(r, "Apart_on") == 0, "Apart fired without block_same_complex_binding");
      double const bound = final_value(r, "Bound_on");
      check(std::abs(bound - expected) < tol,
            "Bound (control) should fire without block_same_complex_binding: " +
                std::to_string(bound));
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: exception: %s\n", e.what());
    ++g_failures;
  }

  if (g_failures == 0)
    std::fprintf(stdout, "products_never_form_test: OK\n");
  return g_failures == 0 ? 0 : 1;
}
