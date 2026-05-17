// Integration test for species enumeration + `.species` output
// (issue #9 §2, plan §6 step 4) and the cached-incremental label layer
// (plan §6 step 5, decision #6).  Exercises the full batch pipeline:
// AgentPool -> extract_complex -> canonical_label -> isomorphism dedup
// -> sorted SpeciesRow list -> BNG-format `.species` file.
//
// Run against three models:
//   A_plus_A              -- A(a)+A(a)<->A(a!1).A(a!1): monomer plus a
//                            symmetric homodimer, so the pool walk hits
//                            the canonical individualization search.
//                            Merge-heavy (and unbind/split).
//   ss_symmetric_homopoly -- P(s,s) self-binding: larger symmetric
//                            chains/rings stress the canonicalizer;
//                            merge- and split-heavy.
//   ft_ring_closure       -- A/B binding with ring closure + reverse
//                            rules: exercises cycle-bond add_bond,
//                            split-on-unbind, and ring formation.
//
// argv: <A_plus_A.xml> <ss_symmetric_homopoly.xml> <ft_ring_closure.xml>
//
// In Debug/ASan builds, enumerate_species cross-checks each complex's
// cached canonical label against a from-scratch recompute and aborts on
// mismatch (the decision-#6 invariant).  test_cached_label_invariant
// calls enumerate_species repeatedly across many simulate() segments so
// that cross-check runs against a cache that has been populated,
// invalidated by merge/split/state events, and lazily recomputed
// mid-run — a missed invalidation hook surfaces there as an abort.

#include "rulemonkey/simulator.hpp"

#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

using rulemonkey::RuleMonkeySimulator;
using rulemonkey::SpeciesRow;

namespace {

int g_failures = 0;

void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

// Number of molecules in a canonical species string = the count of
// `.`-separated molecule blocks.  BNGL identifiers, states and bond
// labels never contain `.`, so it only ever separates molecules.
int molecule_count_of(const std::string& s) {
  int n = 1;
  for (char const c : s)
    if (c == '.')
      ++n;
  return n;
}

// Right after initialize(), before any SSA event, the census must be
// exactly the seed species — a deterministic check independent of seed.
void test_seed_species_exact(const std::string& aa_xml) {
  RuleMonkeySimulator sim(aa_xml);
  sim.initialize(42);
  auto rows = sim.enumerate_species();
  check(rows.size() == 1, "A_plus_A seed census has exactly one species");
  if (rows.size() == 1) {
    check(rows[0].species == "A(a)", "seed species is A(a), got '" + rows[0].species + "'");
    check(rows[0].count == 1000, "seed species count is A_tot = 1000");
  }
}

// After a run: rows sorted and distinct, counts positive, molecules
// conserved, and every dimer collapsed onto one canonical row.
void test_run_conservation_and_grouping(const std::string& aa_xml) {
  RuleMonkeySimulator sim(aa_xml);
  sim.initialize(42);
  sim.simulate(0.0, 10.0, 10);
  auto rows = sim.enumerate_species();

  long total_mol = 0;
  for (size_t i = 0; i < rows.size(); ++i) {
    check(rows[i].count > 0, "species count is positive");
    check(!rows[i].species.empty(), "species string is non-empty");
    if (i)
      check(rows[i - 1].species < rows[i].species, "rows are sorted and distinct by species");
    total_mol += rows[i].count * molecule_count_of(rows[i].species);
    // A_plus_A can only form monomers and the symmetric homodimer; all
    // dimers — however the two molecules were ordered in the pool —
    // must canonicalize to the one string.
    check(rows[i].species == "A(a)" || rows[i].species == "A(a!1).A(a!1)",
          "A_plus_A species is the monomer or the canonical homodimer, got '" + rows[i].species +
              "'");
  }
  check(total_mol == 1000, "molecule conservation: sum of count*size == A_tot (1000)");
}

// write_species_file round-trips: the file's data lines must reproduce
// enumerate_species() exactly (sorted, two-space separator, integer
// count), under a `#` comment header readable by BNG2.pl readNFspecies.
void test_species_file_roundtrip(const std::string& aa_xml) {
  RuleMonkeySimulator sim(aa_xml);
  sim.initialize(42);
  sim.simulate(0.0, 10.0, 10);
  auto rows = sim.enumerate_species();

  const std::string path = "species_enumeration_test_scratch.species";
  sim.write_species_file(path);

  std::ifstream in(path);
  check(in.is_open(), "written species file can be reopened");
  std::vector<SpeciesRow> parsed;
  int header_lines = 0;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    if (line[0] == '#') {
      ++header_lines;
      continue;
    }
    // `<pattern>  <count>` — the pattern has no spaces, so the first
    // space starts the separator gap.
    size_t const sp = line.find(' ');
    check(sp != std::string::npos, "data line has a pattern/count separator");
    if (sp == std::string::npos)
      continue;
    parsed.push_back(SpeciesRow{line.substr(0, sp), std::stol(line.substr(sp))});
  }
  in.close();
  std::remove(path.c_str());

  check(header_lines >= 1, "species file carries a comment header");
  check(parsed.size() == rows.size(), "file row count matches enumerate_species()");
  for (size_t i = 0; i < parsed.size() && i < rows.size(); ++i) {
    check(parsed[i].species == rows[i].species,
          "round-trip species string at row " + std::to_string(i));
    check(parsed[i].count == rows[i].count, "round-trip count at row " + std::to_string(i));
  }
}

// enumerate_species / write_species_file require a live session.
void test_no_session_throws(const std::string& aa_xml) {
  RuleMonkeySimulator sim(aa_xml);
  bool threw = false;
  try {
    (void)sim.enumerate_species();
  } catch (const std::exception&) {
    threw = true;
  }
  check(threw, "enumerate_species() without a session throws");

  threw = false;
  try {
    sim.write_species_file("should_not_be_created.species");
  } catch (const std::exception&) {
    threw = true;
  }
  check(threw, "write_species_file() without a session throws");
}

// A symmetric multi-molecule model: P(s,s) self-binding builds long
// symmetric chains and rings.  The pipeline must canonicalize them all,
// conserve molecules, and keep rows sorted/distinct.
void test_symmetric_pipeline(const std::string& hp_xml) {
  RuleMonkeySimulator sim(hp_xml);
  sim.initialize(42);
  sim.simulate(0.0, 50.0, 10);
  auto rows = sim.enumerate_species();

  check(!rows.empty(), "homopolymer census is non-empty");
  long total_mol = 0;
  bool has_multi = false;
  for (size_t i = 0; i < rows.size(); ++i) {
    check(rows[i].count > 0, "homopolymer species count is positive");
    check(!rows[i].species.empty(), "homopolymer species string is non-empty");
    if (i)
      check(rows[i - 1].species < rows[i].species, "homopolymer rows are sorted and distinct");
    int const n = molecule_count_of(rows[i].species);
    total_mol += rows[i].count * n;
    if (n > 1)
      has_multi = true;
  }
  check(total_mol == 200, "homopolymer conservation: sum of count*size == P_tot (200)");
  check(has_multi, "homopolymer run formed at least one multi-molecule species");
}

// Cached-incremental label invariant (plan decision #6, step 5).  Drive
// many simulate() segments, calling enumerate_species between each.  In
// Debug/ASan builds every enumerate_species runs the decision-#6
// cross-check (cached label == from-scratch recompute) over the current
// cache: complexes untouched since the previous segment are validated
// against their *stale* cache entry, so a structural mutator that fails
// to dirty an edited complex aborts the process here.  In Release the
// check compiles out and this still exercises repeated batch sweeps
// interleaved with simulation.  Independent of the cross-check, the test
// asserts molecule conservation and sorted/distinct rows every segment.
void test_cached_label_invariant(const std::string& xml, const std::string& tag, double seg_len,
                                 long expected_mol_total) {
  RuleMonkeySimulator sim(xml);
  sim.initialize(7);

  double t = 0.0;
  for (int seg = 0; seg < 24; ++seg) {
    auto rows = sim.enumerate_species();
    long total = 0;
    for (size_t i = 0; i < rows.size(); ++i) {
      check(rows[i].count > 0, tag + ": species count positive");
      total += rows[i].count * molecule_count_of(rows[i].species);
      if (i)
        check(rows[i - 1].species < rows[i].species, tag + ": rows sorted and distinct");
    }
    check(total == expected_mol_total,
          tag + ": molecule conservation at segment " + std::to_string(seg) + " (got " +
              std::to_string(total) + ", expected " + std::to_string(expected_mol_total) + ")");
    double const next = t + seg_len;
    sim.simulate(t, next, 1);
    t = next;
  }
  // A final sweep after the last simulate segment so the cross-check
  // also validates the cache state left by the trailing events.
  auto rows = sim.enumerate_species();
  long total = 0;
  for (const auto& r : rows)
    total += r.count * molecule_count_of(r.species);
  check(total == expected_mol_total, tag + ": molecule conservation after final segment");
}

// load_state must clear the label cache wholesale (the restored pool has
// all-new complex ids).  Save mid-run, load into a fresh simulator, then
// enumerate — in Debug/ASan the cross-check would abort if a stale entry
// from before the load survived.
void test_cache_survives_save_load(const std::string& xml, const std::string& tag) {
  const std::string path = "species_enumeration_test_state.scratch";
  RuleMonkeySimulator sim(xml);
  sim.initialize(7);
  sim.simulate(0.0, 5.0, 5);
  (void)sim.enumerate_species(); // populate the label cache
  sim.save_state(path);

  RuleMonkeySimulator loaded(xml);
  loaded.initialize(99);
  loaded.simulate(0.0, 2.0, 2); // give the fresh sim its own cache state
  (void)loaded.enumerate_species();
  loaded.load_state(path); // must drop that cache wholesale
  auto rows = loaded.enumerate_species();
  std::remove(path.c_str());
  check(!rows.empty(), tag + ": species census non-empty after load_state");
  // The loaded session resumes at the save point (t=5.0); keep editing
  // the pool post-load so the cross-check also covers the rebuilt cache.
  loaded.simulate(5.0, 9.0, 4);
  (void)loaded.enumerate_species();
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::fprintf(stderr, "usage: species_enumeration_test <A_plus_A.xml> "
                         "<ss_symmetric_homopoly.xml> <ft_ring_closure.xml>\n");
    return 2;
  }
  const std::string aa_xml = argv[1];
  const std::string hp_xml = argv[2];
  const std::string ring_xml = argv[3];

  try {
    test_seed_species_exact(aa_xml);
    test_run_conservation_and_grouping(aa_xml);
    test_species_file_roundtrip(aa_xml);
    test_no_session_throws(aa_xml);
    test_symmetric_pipeline(hp_xml);
    // Cached-incremental layer: merge-heavy, split-heavy, and ring-forming
    // models, each swept many times mid-run.  Seed totals: A_plus_A
    // A_tot=1000; ss_symmetric_homopoly P_tot=200; ft_ring_closure
    // A_tot=B_tot=60 -> 120 molecules.
    test_cached_label_invariant(aa_xml, "A_plus_A", 2.0, 1000);
    test_cached_label_invariant(hp_xml, "ss_symmetric_homopoly", 5.0, 200);
    test_cached_label_invariant(ring_xml, "ft_ring_closure", 5.0, 120);
    test_cache_survives_save_load(hp_xml, "ss_symmetric_homopoly");
    test_cache_survives_save_load(ring_xml, "ft_ring_closure");
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: species enumeration + .species output all pass\n");
  return 0;
}
