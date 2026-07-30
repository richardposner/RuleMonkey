// Corpus-wide guard on the symbolic-baseline gate (issue #23).
//
// The set_param cascade re-derives parameters from `<Parameter expr=>`,
// because BNG2's `<Parameter value=>` is already collapsed to a number and
// therefore cannot propagate an override.  But `value` is what RuleMonkey
// resolves at load time — deliberately, for NFsim parity — and BNG2 does
// not always write the two attributes to the same precision:
//
//   <Parameter id="NA" value="6.0221408e+23" expr="6.02214076e23"/>
//
// So the cascade only accepts an expr-derived value where it actually
// MOVES the parameter off its no-override baseline; anything the override
// does not reach keeps its loaded `value`.  Without that gate, a single
// set_param anywhere in a model would silently re-round every derived
// quantity in it — including every bimolecular rate constant that divides
// by Avogadro's number.
//
// This test asserts the gate holds across every model XML in the tree.
// For each model it issues a NO-OP override (set_param on the first
// declared parameter, to the value that parameter already has), which is
// enough to switch the whole cascade onto the `expr` source, and then
// requires every parameter to be bit-identical to the un-overridden load.
// A single model whose `value=` and `expr=` disagree beyond the gate's
// reach would fail here rather than silently shifting a trajectory.
//
// Load-only, so the whole corpus sweeps in about a second.

#include "rulemonkey/simulator.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <tests-dir>\n", argv[0]);
    return 2;
  }

  // Deterministic order: directory iteration order is unspecified.
  std::vector<std::string> xmls;
  std::error_code ec;
  for (const auto& entry : std::filesystem::recursive_directory_iterator(argv[1], ec)) {
    if (entry.is_regular_file() && entry.path().extension() == ".xml")
      xmls.push_back(entry.path().string());
  }
  if (ec) {
    std::fprintf(stderr, "ERROR: cannot walk '%s': %s\n", argv[1], ec.message().c_str());
    return 2;
  }
  std::sort(xmls.begin(), xmls.end());
  if (xmls.size() < 50) {
    std::fprintf(stderr,
                 "ERROR: only %zu XML models found under '%s' — expected the full "
                 "corpus; the test would pass vacuously\n",
                 xmls.size(), argv[1]);
    return 2;
  }

  int swept = 0;
  for (const auto& xml : xmls) {
    try {
      const rulemonkey::RuleMonkeySimulator base(xml);
      rulemonkey::RuleMonkeySimulator over(xml);
      const auto names = base.parameter_names();
      if (names.empty())
        continue;
      ++swept;

      over.set_param(names.front(), base.get_parameter(names.front()));

      for (const auto& name : names) {
        const double before = base.get_parameter(name);
        const double after = over.get_parameter(name);
        if (before != after) {
          std::fprintf(stderr,
                       "FAIL: %s: a no-op override moved parameter '%s' from %.17g to %.17g "
                       "(value=/expr= drift escaping the symbolic-baseline gate)\n",
                       xml.c_str(), name.c_str(), before, after);
          ++g_failures;
        }
      }
    } catch (const std::exception& e) {
      std::fprintf(stderr, "FAIL: %s: %s\n", xml.c_str(), e.what());
      ++g_failures;
    }
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed across %d models\n", g_failures, swept);
    return 1;
  }
  std::fprintf(stderr,
               "OK: %d models swept, no parameter drifts when the set_param cascade "
               "switches to the expr= source\n",
               swept);
  return 0;
}
