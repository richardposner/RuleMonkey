// What -no-bscb still lets a bimolecular rule do inside one complex, once the
// `+` between product patterns is tested whatever block_same_complex_binding
// says (issue #123).
//
// The BNG2-oracle model tests/bng_oracle/models/plus_no_bscb.bngl pins the
// other half under -no-bscb: a rule whose `+` cannot hold never fires. BNG2
// cannot referee this half, since its network never lets one complex supply
// both reactant patterns, so it is pinned here. The fixture,
// plus_without_bscb_model.xml, holds one E.S complex and one F.G complex and
// nothing else, so every pair either rule can draw lies inside one complex.
//
//   * `Mod` (`E(s) + S(p~0) -> E(s) + S(p~1)`) deletes nothing. Under
//     -no-bscb it fires on the E and S of one complex: the `+` asks what the
//     rule parts, and it parts nothing. This is NFsim's reading.
//   * `Close` (`F(x) + G(y) -> F(x!1).G(y!1)`) closes a ring, the job
//     -no-bscb exists for.
//
// Under the default, neither fires. k * t_end = 20, so a rule that can fire
// has fired but for a chance of e^-20.

#include "rulemonkey/simulator.hpp"

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
    std::fprintf(stderr, "Usage: plus_without_bscb_test <plus_without_bscb_model.xml>\n");
    return 2;
  }

  try {
    for (bool const bscb : {true, false}) {
      std::string const mode =
          bscb ? "under block_same_complex_binding" : "without block_same_complex_binding";
      for (std::uint64_t const seed : {1, 2, 3}) {
        rulemonkey::RuleMonkeySimulator sim(argv[1]);
        sim.set_block_same_complex_binding(bscb);
        auto r = sim.run(rulemonkey::TimeSpec{0.0, 5.0, 2}, seed);
        double const want = bscb ? 0.0 : 1.0;
        check(final_value(r, "Mod_on") == want,
              "Mod " + mode + ", seed " + std::to_string(seed) + ": expected " +
                  std::to_string(want) + ", observed " + std::to_string(final_value(r, "Mod_on")));
        check(final_value(r, "Close_on") == want, "Close " + mode + ", seed " +
                                                      std::to_string(seed) + ": expected " +
                                                      std::to_string(want) + ", observed " +
                                                      std::to_string(final_value(r, "Close_on")));
      }
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: exception: %s\n", e.what());
    ++g_failures;
  }

  if (g_failures == 0)
    std::fprintf(stdout, "plus_without_bscb_test: OK\n");
  return g_failures == 0 ? 0 : 1;
}
