// Error-path regression test.  Pins down that the public API throws
// std::runtime_error (not std::exception, not silent failure, not
// crash) for the diagnostic surfaces it documents:
//
//   1. Constructor with a path that does not exist.
//   2. Constructor with malformed XML on disk.
//   3. set_param with an undeclared parameter name.
//   4. Mutators called while a session is active:
//        set_param, clear_param_overrides, set_molecule_limit,
//        set_block_same_complex_binding.
//
// These paths exist in simulator.cpp/engine.cpp but were previously
// only exercised indirectly by the corpus parity tests.  A regression
// that swallowed an XML-load error or accepted a stale set_param
// call after run() would never have surfaced as a parity failure;
// here it surfaces as a one-line diff.

#include "rulemonkey/simulator.hpp"

#include <cstdio>
#include <fstream>
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

// True iff calling `fn()` throws std::runtime_error.  Catches
// std::exception too, but reports the type so a regression that
// downgrades runtime_error to a plain exception is visible.
template <typename F> bool throws_runtime_error(F&& fn, const std::string& tag) {
  try {
    fn();
  } catch (const std::runtime_error&) {
    return true;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "  %s: caught std::exception (not runtime_error): %s\n", tag.c_str(),
                 e.what());
    return false;
  } catch (...) {
    std::fprintf(stderr, "  %s: caught non-std exception\n", tag.c_str());
    return false;
  }
  std::fprintf(stderr, "  %s: did not throw\n", tag.c_str());
  return false;
}

void test_missing_file() {
  // Path that cannot exist.  The simulator should refuse to construct.
  check(throws_runtime_error(
            []() { rulemonkey::RuleMonkeySimulator sim("/nonexistent/path/__no_such_xml__.xml"); },
            "missing_file"),
        "constructor must throw runtime_error on missing file");
}

void test_malformed_xml(const std::string& tmp_path) {
  // Drop a file that is syntactically not XML.  Constructor must fail
  // loudly rather than producing an empty / partially-loaded model.
  {
    std::ofstream f(tmp_path);
    f << "this is not xml { not even close }\n";
  }
  check(throws_runtime_error([&]() { rulemonkey::RuleMonkeySimulator sim(tmp_path); },
                             "malformed_xml"),
        "constructor must throw runtime_error on malformed XML");
  std::remove(tmp_path.c_str());
}

void test_set_param_unknown(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  check(
      throws_runtime_error([&]() { sim.set_param("__no_such_param__", 1.0); }, "set_param_unknown"),
      "set_param must throw on undeclared parameter name");
}

void test_mutators_during_session(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(/*seed=*/1);
  check(sim.has_session(), "session should be active after initialize()");

  check(throws_runtime_error([&]() { sim.set_param("kp", 1.0); }, "set_param_active"),
        "set_param must throw while a session is active");

  check(
      throws_runtime_error([&]() { sim.clear_param_overrides(); }, "clear_param_overrides_active"),
      "clear_param_overrides must throw while a session is active");

  check(throws_runtime_error([&]() { sim.set_molecule_limit(100); }, "set_molecule_limit_active"),
        "set_molecule_limit must throw while a session is active");

  check(throws_runtime_error([&]() { sim.set_block_same_complex_binding(false); },
                             "set_block_same_complex_binding_active"),
        "set_block_same_complex_binding must throw while a session is active");

  // After destroy_session, the same calls should succeed.
  sim.destroy_session();
  check(!sim.has_session(), "session should be inactive after destroy_session()");
  try {
    sim.set_param("kp", 1.0);
    sim.clear_param_overrides();
    sim.set_molecule_limit(1000);
    sim.set_block_same_complex_binding(true);
  } catch (const std::exception& e) {
    check(false, std::string("post-destroy mutators must not throw: ") + e.what());
  }
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <A_plus_A.xml> <tmp_writable_path.xml>\n", argv[0]);
    return 2;
  }
  const std::string xml = argv[1];
  const std::string tmp_path = argv[2];

  test_missing_file();
  test_malformed_xml(tmp_path);
  test_set_param_unknown(xml);
  test_mutators_during_session(xml);

  if (g_failures != 0) {
    std::fprintf(stderr, "error_paths_test: %d failure(s)\n", g_failures);
    return 1;
  }
  std::printf("error_paths_test: OK\n");
  return 0;
}
