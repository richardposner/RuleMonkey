// Regression test: save_state / load_state round-trip + schema-mismatch refusal.
//
// Two contracts:
//   1. round-trip: save mid-run, load into a fresh simulator, the
//      continuation trajectory matches what an uninterrupted single
//      run produces from the same seed (deterministic by construction).
//   2. fingerprint guard: a state file saved against XML A must NOT
//      load into a simulator constructed from a structurally different
//      XML B.  Pre-fix the load was silently accepted and the loaded
//      pool's component indices referred to A's schema, producing
//      corrupt observables on the next run.

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

// Round-trip: a save-at-t/load continuation reproduces an uninterrupted run.
//
// The deterministic equality follows from save_state capturing the
// RNG state immediately before the next dt draw and load_state
// re-seeding from that snapshot — so the second leg redraws bit-
// identical numbers.
void test_round_trip(const std::string& xml) {
  const std::uint64_t seed = 7;
  const double t_mid = 5.0;
  const double t_end = 10.0;
  const int n_pts = 11;

  // Reference: one uninterrupted run from 0..t_end.
  rulemonkey::RuleMonkeySimulator ref(xml);
  ref.set_block_same_complex_binding(true);
  rulemonkey::Result const ref_full = ref.run({0.0, t_end, n_pts}, seed);
  double const ref_aa_end = final_obs(ref_full, "AA_1");

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
  rulemonkey::Result const tail = s2.simulate(t_mid, t_end, 6);
  std::filesystem::remove(tmp);

  double const split_aa_end = final_obs(tail, "AA_1");
  check(ref_aa_end == split_aa_end,
        "save/load round-trip: continuation final AA_1 should match uninterrupted run "
        "(ref=" +
            std::to_string(ref_aa_end) + " vs split=" + std::to_string(split_aa_end) + ")");
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
  std::fprintf(stderr, "OK: save/load round-trip preserves trajectory; XML mismatch refused; "
                       "V1 format rejected\n");
  return 0;
}
