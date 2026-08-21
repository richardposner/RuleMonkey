// The per-rule side-table accounting reported by `[RM tables]` (issue #71).
//
// #71 is about what a rule's per-molecule tables COST: they are sized by the
// molecule arena rather than by the population the rule can see, so a model
// holds a (rule count x arena x row width) product in them.  #72 attacked
// the row width and left the sizing, on the grounds that the case for
// changing it is residency on real models — which is a measurement, and the
// engine now reports the number that measurement needs.
//
// This pins that report.  Three fixtures, all of them models already in the
// tree for other tests, chosen so that between them every bucket of the
// accounting is represented:
//
//   pool_churn_model      the narrow row, the wide row, and — the point of
//                         #68 — a rule with no reactant seed molecule
//                         (`0 -> W()`), which holds no table at all and so
//                         must not appear in the report.
//   trimolecular_model    an n-ary rule, whose per-slot count tables are a
//                         different set of vectors from `mol_data`, and
//                         whose 400-molecule seed type is over
//                         FENWICK_THRESHOLD, so its samplers are allocated.
//   rule_table_shape_model  #72's own fixture: one rule per wide-half group
//                         (shared-component split, local rate law,
//                         pure-context slot) plus the shapes with no wide
//                         half.
//
// What the assertions pin, and what they do not.  Each rule's reported
// per-molecule bytes must be an exact multiple of its reported row count,
// and the multiplier must be one of the two widths a two-slot rule can
// have: `kNarrowRow` for a rule that reads only `PerMolRuleData`, and
// `kNarrowRow + kWideRow` for one that also carries a `PerMolRuleAux`.  A
// failure here is not necessarily a bug — it is what a field added to
// either struct looks like — but it is the one change #71 says has to be
// made deliberately, since a byte on the row is a byte per molecule per
// rule of resident memory for the whole run.  Update the constant and say
// what the field bought.
//
// What no assertion here can catch is a NEW arena-sized per-rule table
// added and not added to `rule_table_bytes`: it would go unreported and
// every width below would still check out.  The report is an accounting of
// what the engine knows it holds, not a measurement of the heap.
//
// The report is stderr, gated on RM_PRINT_TABLES, so the fixture sets the
// variable before the first run — the engine reads it once into a
// function-local static — and captures stderr the way
// negative_rate_clamp_test does.

#include "rulemonkey/simulator.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::fprintf(stdout, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

// sizeof(PerMolRuleData) and sizeof(PerMolRuleAux) as of #72.  Deliberately
// spelled out rather than derived: the engine's structs are private to
// engine.cpp, and a test that recomputed the width from the same source it
// is checking would pin nothing.
constexpr unsigned long long kNarrowRow = 12; // two int counts + the P1 cache flag
constexpr unsigned long long kWideRow = 64;   // seven doubles + two complex ids

// One `[RM tables]` block: the summary line plus a line per rule holding a
// table.
struct Block {
  unsigned long long arena = 0;
  unsigned long long rules = 0;
  unsigned long long tabled = 0;
  unsigned long long bytes = 0;
  unsigned long long mol = 0;
  unsigned long long sampler = 0;
  unsigned long long reach_bytes = 0;
  struct Rule {
    std::string id;
    unsigned long long rows = 0;
    unsigned long long reach = 0;
    unsigned long long mol = 0;
    unsigned long long sampler = 0;
  };
  std::vector<Rule> per_rule;
};

std::vector<Block> parse_blocks(const std::string& blob) {
  std::vector<Block> blocks;
  std::istringstream lines(blob);
  std::string line;
  while (std::getline(lines, line)) {
    Block b;
    double mb_unused = 0;
    // NOLINTNEXTLINE(cert-err34-c) — fixture parsing of our own output.
    if (std::sscanf(line.c_str(),
                    "[RM tables] arena=%llu rules=%llu tabled=%llu bytes=%llu (%lf MB) mol=%llu "
                    "sampler=%llu reach_bytes=%llu",
                    &b.arena, &b.rules, &b.tabled, &b.bytes, &mb_unused, &b.mol, &b.sampler,
                    &b.reach_bytes) != 8)
      continue;
    // Per-rule detail runs until the next line that is not one.
    std::streampos next = lines.tellg();
    while (std::getline(lines, line)) {
      Block::Rule r;
      char id[256] = {0};
      // NOLINTNEXTLINE(cert-err34-c)
      if (std::sscanf(line.c_str(),
                      "  %255[^ ] (%*[^)]): rows=%llu reach=%llu mol=%llu sampler=%llu",
                      static_cast<char*>(id), &r.rows, &r.reach, &r.mol, &r.sampler) != 5) {
        lines.clear();
        lines.seekg(next);
        break;
      }
      r.id = id;
      b.per_rule.push_back(r);
      next = lines.tellg();
    }
    blocks.push_back(b);
  }
  return blocks;
}

// The summary has to be the sum of the detail, and `tabled` has to be how
// many detail lines there are: a sweep reads the summary and trusts it.
void check_self_consistent(const Block& b, const std::string& what) {
  unsigned long long mol = 0;
  unsigned long long sampler = 0;
  for (const auto& r : b.per_rule) {
    mol += r.mol;
    sampler += r.sampler;
  }
  check(b.tabled == b.per_rule.size(), what + ": tabled=" + std::to_string(b.tabled) + " but " +
                                           std::to_string(b.per_rule.size()) + " detail lines");
  check(b.mol == mol, what + ": summary mol=" + std::to_string(b.mol) + " but detail sums to " +
                          std::to_string(mol));
  check(b.sampler == sampler, what + ": summary sampler=" + std::to_string(b.sampler) +
                                  " but detail sums to " + std::to_string(sampler));
  check(b.bytes == b.mol + b.sampler, what + ": bytes != mol + sampler");
  // Sizing by reach can only ever shorten a table, never lengthen one.
  check(b.reach_bytes <= b.bytes, what + ": reach_bytes=" + std::to_string(b.reach_bytes) +
                                      " exceeds bytes=" + std::to_string(b.bytes));
}

// Every two-slot rule's row is one of the two widths, and nothing is
// charged for a fraction of a row.  Returns how many carried the wide half.
int check_two_slot_widths(const Block& b, const std::string& what) {
  int wide = 0;
  for (const auto& r : b.per_rule) {
    check(r.rows > 0, what + "/" + r.id + ": table with no rows");
    if (r.rows == 0)
      continue;
    unsigned long long const width = r.mol / r.rows;
    check(width * r.rows == r.mol, what + "/" + r.id + ": mol=" + std::to_string(r.mol) +
                                       " is not a whole multiple of " + std::to_string(r.rows) +
                                       " rows");
    check(width == kNarrowRow || width == kNarrowRow + kWideRow,
          what + "/" + r.id + ": row width " + std::to_string(width) + " is neither " +
              std::to_string(kNarrowRow) + " (narrow) nor " +
              std::to_string(kNarrowRow + kWideRow) + " (narrow + wide) — a field was added to " +
              "PerMolRuleData or PerMolRuleAux, which costs that many bytes per molecule per " +
              "rule of resident memory (issue #71)");
    if (width == kNarrowRow + kWideRow)
      ++wide;
  }
  return wide;
}

rulemonkey::Result run(const std::string& xml, double t_end) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  rulemonkey::TimeSpec ts;
  ts.t_start = 0.0;
  ts.t_end = t_end;
  ts.n_points = 10;
  return sim.run(ts, /*seed=*/1);
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 5) {
    std::fprintf(stdout, "Usage: rule_table_footprint_test <pool_churn.xml> <trimolecular.xml>"
                         " <rule_table_shape.xml> <stderr_capture_path>\n");
    return 2;
  }
  const std::string churn_xml = argv[1];
  const std::string trimol_xml = argv[2];
  const std::string shape_xml = argv[3];
  const std::string err_path = argv[4];

  // Before the first run: the engine latches the gate into a function-local
  // static the first time a run finishes.  MSVC has no `setenv`, and the
  // Windows CI leg builds and runs this suite.
#ifdef _MSC_VER
  _putenv_s("RM_PRINT_TABLES", "1");
#else
  // NOLINTNEXTLINE(concurrency-mt-unsafe) — single-threaded fixture setup.
  setenv("RM_PRINT_TABLES", "1", 1);
#endif

  if (std::freopen(err_path.c_str(), "w", stderr) == nullptr) {
    std::fprintf(stdout, "FAIL: could not redirect stderr to %s\n", err_path.c_str());
    return 1;
  }
  try {
    run(churn_xml, 1.0);
    run(trimol_xml, 1.0);
    run(shape_xml, 1.0);
  } catch (const std::exception& e) {
    std::fflush(stderr);
    std::fprintf(stdout, "FAIL: sim.run threw: %s\n", e.what());
    return 1;
  }
  std::fflush(stderr);

  std::ifstream const ef(err_path);
  std::stringstream ss;
  ss << ef.rdbuf();
  const std::string blob = ss.str();

  auto const blocks = parse_blocks(blob);
  if (blocks.size() != 3) {
    std::fprintf(stdout, "--- captured stderr ---\n%s--- end ---\n", blob.c_str());
    std::fprintf(stdout, "FAIL: expected one [RM tables] block per run, got %zu\n", blocks.size());
    return 1;
  }

  // ---- pool_churn_model -------------------------------------------------
  // Five rules: X->Y, Y->X, the Z homodimer pair, and `0 -> W()`.  The last
  // has no reactant seed molecule, so #68 leaves it no table and it must not
  // be in the report at all — the accounting has to show that saving, since
  // it is the one #71's direction 3 is the general case of.
  {
    Block const& b = blocks[0];
    check_self_consistent(b, "pool_churn");
    check(b.rules == 5, "pool_churn: expected 5 rules, got " + std::to_string(b.rules));
    check(b.tabled == 4,
          "pool_churn: expected 4 rules holding a table (`0 -> W()` holds none, #68), got " +
              std::to_string(b.tabled));
    int const wide = check_two_slot_widths(b, "pool_churn");
    check(wide == 1, "pool_churn: exactly one rule (the Z homodimer forward) should carry the "
                     "wide half, got " +
                         std::to_string(wide));
    check(b.sampler == 0,
          "pool_churn: every type is under FENWICK_THRESHOLD, so no sampler should be allocated");
    for (const auto& r : b.per_rule)
      check(r.reach > 0, "pool_churn/" + r.id + ": a rule holding a table reaches no molecule");
  }

  // ---- trimolecular_model -----------------------------------------------
  // One n-ary rule, `A(s) + A(s) + A(s) -> P()`, seeded with 400 A.  Its
  // per-molecule state is three `counts` vectors of int rather than a
  // `mol_data`, and 400 is over FENWICK_THRESHOLD (300), so all three slots
  // get a sampler — three trees of (arena at session build + 1) doubles,
  // the Fenwick's 1-indexed layout.
  {
    Block const& b = blocks[1];
    check_self_consistent(b, "trimolecular");
    check(b.rules == 1, "trimolecular: expected 1 rule, got " + std::to_string(b.rules));
    check(b.tabled == 1, "trimolecular: expected 1 rule holding a table");
    if (b.per_rule.size() == 1) {
      Block::Rule const& r = b.per_rule[0];
      check(r.mol == r.rows * 3 * sizeof(int), "trimolecular: expected three int slot tables of " +
                                                   std::to_string(r.rows) +
                                                   " rows, got mol=" + std::to_string(r.mol));
      unsigned long long const expect_sampler = 3ULL * (400 + 1) * sizeof(double);
      check(r.sampler == expect_sampler,
            "trimolecular: expected three Fenwick trees over the 400-molecule seed arena (" +
                std::to_string(expect_sampler) + " bytes), got " + std::to_string(r.sampler));
      // The 400 A are seeded first, so they hold ids 0..399 and every P the
      // rule makes lands above them: the rule reaches 400 of an arena that
      // its own products push past that.  The model is small enough to pin
      // that exactly, which is what makes `reach` a checked number rather
      // than a reported one.
      check(r.reach == 400,
            "trimolecular: expected reach 400 (the seeded A), got " + std::to_string(r.reach));
    }
  }

  // ---- rule_table_shape_model -------------------------------------------
  // #72's fixture, whose nine expanded rules cover all three wide-half
  // groups.  Exactly three carry one: the A homodimer forward (the
  // shared-component split), `Lrule` (a local rate law), and the E+S
  // catalysis (a pure-context slot).  The reverse of the homodimer is
  // unimolecular and the four makers and the D decay are plain, so those six
  // stay narrow.
  {
    Block const& b = blocks[2];
    check_self_consistent(b, "rule_table_shape");
    check(b.rules == 9, "rule_table_shape: expected 9 rules, got " + std::to_string(b.rules));
    check(b.tabled == 9, "rule_table_shape: every rule seeds on a molecule, so all 9 hold a table");
    int const wide = check_two_slot_widths(b, "rule_table_shape");
    check(wide == 3, "rule_table_shape: exactly three shapes read the wide half (homodimer "
                     "forward, local rate law, pure-context slot), got " +
                         std::to_string(wide));
  }

  if (g_failures > 0)
    std::fprintf(stdout, "--- captured stderr ---\n%s--- end ---\n", blob.c_str());
  std::fprintf(stdout, "%s (%d failures)\n", g_failures == 0 ? "PASS" : "FAILED", g_failures);
  return g_failures == 0 ? 0 : 1;
}
