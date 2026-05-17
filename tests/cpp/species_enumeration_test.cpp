// Integration test for species enumeration + `.species` output
// (issue #9 §2, plan §6 step 4).  Exercises the full batch pipeline:
// AgentPool -> extract_complex -> canonical_label -> isomorphism dedup
// -> sorted SpeciesRow list -> BNG-format `.species` file.
//
// Run against two models:
//   A_plus_A              -- A(a)+A(a)<->A(a!1).A(a!1): monomer plus a
//                            symmetric homodimer, so the pool walk hits
//                            the canonical individualization search.
//   ss_symmetric_homopoly -- P(s,s) self-binding: larger symmetric
//                            chains/rings stress the canonicalizer.
//
// argv: <A_plus_A.xml> <ss_symmetric_homopoly.xml>

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

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stderr,
                 "usage: species_enumeration_test <A_plus_A.xml> <ss_symmetric_homopoly.xml>\n");
    return 2;
  }
  const std::string aa_xml = argv[1];
  const std::string hp_xml = argv[2];

  try {
    test_seed_species_exact(aa_xml);
    test_run_conservation_and_grouping(aa_xml);
    test_species_file_roundtrip(aa_xml);
    test_no_session_throws(aa_xml);
    test_symmetric_pipeline(hp_xml);
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
