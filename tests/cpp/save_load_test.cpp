// Regression test: save_state / load_state round-trip + schema-mismatch refusal.
//
// Three contracts:
//   1. round-trip: a session saved mid-run and loaded into a fresh
//      simulator continues exactly as the same session continued in
//      process: the same samples and the same event count, from the same
//      seed.  It does not reproduce an uninterrupted run.  A call ends on
//      a waiting-time draw that lands past its t_end; that draw is spent,
//      and the next call takes a fresh one (issue #89), whereas a run
//      with no boundary there fires where the draw lands.  Both are
//      exact, and they part ways from the boundary on.  This used to
//      compare against the uninterrupted run, which held only because the
//      RNG was rolled back to before the spent draw so that the next call
//      redrew it, and that reuse undercounted.  What save/load must keep
//      is agreement with the in-process continuation: save_state writes
//      the RNG past the spent draw, so both take the same next number.
//   2. fingerprint guard: a state file saved against XML A must NOT
//      load into a simulator constructed from a structurally different
//      XML B.  Pre-fix the load was silently accepted and the loaded
//      pool's component indices referred to A's schema, producing
//      corrupt observables on the next run.
//   3. the pre-V2 state-file marker is refused outright.

#include "rulemonkey/simulator.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
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

double final_obs(const rulemonkey::Result& r, const std::string& name) {
  for (std::size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == name)
      return r.observable_data[i].back();
  throw std::runtime_error("observable not found: " + name);
}

// Round-trip: a save-at-t/load continuation reproduces the in-process one.
// Compared over the whole tail, not only its last sample.
void test_round_trip(const std::string& xml) {
  const std::uint64_t seed = 7;
  const double t_mid = 5.0;
  const double t_end = 10.0;

  // In process: 0..t_mid, then t_mid..t_end on the same session.
  rulemonkey::RuleMonkeySimulator in_process(xml);
  in_process.set_block_same_complex_binding(true);
  in_process.initialize(seed);
  in_process.simulate(0.0, t_mid, 6);
  rulemonkey::Result const tail_a = in_process.simulate(t_mid, t_end, 6);

  // Split run: 0..t_mid → save → fresh sim → load → t_mid..t_end.
  auto tmp = std::filesystem::temp_directory_path() / "rm_save_load_test.bin";
  {
    rulemonkey::RuleMonkeySimulator s1(xml);
    s1.set_block_same_complex_binding(true);
    s1.initialize(seed);
    s1.simulate(0.0, t_mid, 6);
    s1.save_state(tmp.string());
    s1.destroy_session();
  }

  rulemonkey::RuleMonkeySimulator s2(xml);
  s2.set_block_same_complex_binding(true);
  s2.load_state(tmp.string());
  rulemonkey::Result const tail_b = s2.simulate(t_mid, t_end, 6);
  std::filesystem::remove(tmp);

  check(tail_a.time == tail_b.time, "save/load round-trip: tail sample times differ");
  check(tail_a.observable_data == tail_b.observable_data,
        "save/load round-trip: tail should match the in-process continuation at every sample "
        "(final AA_1 in process=" +
            std::to_string(final_obs(tail_a, "AA_1")) +
            " vs loaded=" + std::to_string(final_obs(tail_b, "AA_1")) + ")");
  // Result::event_count counts the whole session, and the state file
  // carries it, so both tails report the count at t_end.
  check(tail_a.event_count == tail_b.event_count,
        "save/load round-trip: in-process session ended on " + std::to_string(tail_a.event_count) +
            " events, the loaded one on " + std::to_string(tail_b.event_count));
}

// Fingerprint guard: a state saved against XML A must refuse to load
// into a simulator constructed from XML B with a structurally
// different molecule-type schema.  ft_mm_ratelaw introduces molecule
// types absent from A_plus_A, so its fingerprint differs.
void test_fingerprint_mismatch(const std::string& xml_a, const std::string& xml_b) {
  auto tmp = std::filesystem::temp_directory_path() / "rm_save_load_mismatch.bin";

  // Save against XML A.
  {
    rulemonkey::RuleMonkeySimulator sa(xml_a);
    sa.set_block_same_complex_binding(true);
    sa.initialize(/*seed=*/3);
    sa.simulate(0.0, 1.0, 2);
    sa.save_state(tmp.string());
    sa.destroy_session();
  }

  // Load into a simulator constructed from XML B (different schema).
  bool threw = false;
  std::string what;
  try {
    rulemonkey::RuleMonkeySimulator sb(xml_b);
    sb.set_block_same_complex_binding(true);
    sb.load_state(tmp.string());
  } catch (const std::runtime_error& e) {
    threw = true;
    what = e.what();
  }
  std::filesystem::remove(tmp);

  check(threw, "load_state should refuse a state file saved against a different XML "
               "(fingerprint mismatch)");
  if (threw) {
    check(what.find("fingerprint") != std::string::npos,
          "fingerprint-mismatch error message should mention 'fingerprint' (got: " + what + ")");
  }
}

// Reject pre-V2 marker outright.  No production user has shipped V1
// state files (save/load was internal-only through 3.1.x), but we
// document the version transition by failing loudly rather than
// silently accepting an un-fingerprinted file.
void test_v1_marker_rejected(const std::string& xml) {
  auto tmp = std::filesystem::temp_directory_path() / "rm_save_load_v1.bin";
  {
    std::FILE* f = std::fopen(tmp.string().c_str(), "w");
    check(f != nullptr, "could not open temp file for V1 fixture");
    if (f) {
      std::fprintf(f, "RM_STATE_V1\n0.0\n0\n0\n0\n0\nEND\n");
      std::fclose(f);
    }
  }
  bool threw = false;
  std::string what;
  try {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.load_state(tmp.string());
  } catch (const std::runtime_error& e) {
    threw = true;
    what = e.what();
  }
  std::filesystem::remove(tmp);
  check(threw, "load_state should refuse an RM_STATE_V1 file (pre-fingerprint format)");
  if (threw) {
    check(what.find("V1") != std::string::npos,
          "V1-rejection error should mention 'V1' (got: " + what + ")");
  }
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <A_plus_A.xml> <ft_mm_ratelaw.xml>\n", argv[0]);
    return 2;
  }
  const std::string xml_a = argv[1];
  const std::string xml_b = argv[2];

  try {
    test_round_trip(xml_a);
    test_fingerprint_mismatch(xml_a, xml_b);
    test_v1_marker_rejected(xml_a);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: save/load continuation matches in-process; XML mismatch refused; "
                       "V1 format rejected\n");
  return 0;
}
