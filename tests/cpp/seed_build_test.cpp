// Seed-species build (issue #67).
//
// `init_species` used to redo, once per COPY of a seed species, work
// that depends only on the species and the model: the species-XML
// component index -> molecule-type component index mapping (an
// `unordered_map` keyed by component name, built and thrown away per
// molecule), the state-name lookups behind it, the bond endpoints
// resolved through it, and a pair of scratch vectors per copy.  It now
// resolves each seed species once and stamps out the copies from the
// result, and hands the pool the seed totals so its arenas stop growing
// by doubling.  Every copy therefore comes from a template that no
// longer sees the copy it is making, which is exactly the kind of change
// that can leave copy 1 right and copies 2..N wrong.
//
// `seed_build_model.xml` is shaped so each part of that resolution has a
// consequence an exact count can see (see the header comments there and
// in the .bngl).  Every assertion below is a hand-derived integer, not a
// trajectory statistic:
//
//   Pss01 = 3   the three P(s~0,s~1,t~b) monomers: the two `s` slots got
//               DIFFERENT states, so the N-th-occurrence match between
//               the species' duplicate component names and the type's
//               survived.
//   Pss1  = 4   the four dimer P's, both `s` slots in state 1.
//   Ptb   = 3   `t` in state b and unbound — on the monomers only.
//   PtaFree = 0 and PtaBound = Bound = QyBound = 4: the dimer's bond
//               landed on `t`, not on an `s`.  The seed XML lists `t`
//               first, so a bond endpoint resolved without the mapping
//               would land on `s` and flip all four of these.
//   Qx    = 0   the seed species with a count of zero contributes
//               nothing.
//   Pn = 7, Qn = 4, Rn = 0.
//
// The last arm runs the session past the seed.  The model's rules split
// a complex, edit a component state in place, and synthesize a molecule
// from nothing — and that third one is why this file also exists: a
// newborn complex no longer marks itself dirty in the pool's
// canonical-label invalidation set, on the grounds that its id was just
// minted from a counter that never recycles and so cannot be in any
// label cache.  In Debug and ASan builds `enumerate_species` proves
// every cached label against a from-scratch recompute, so this arm is
// where a wrong answer to that would abort.

#include "rulemonkey/simulator.hpp"

#include <cstdio>
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

void check_eq(double got, double want, const std::string& what) {
  check(got == want,
        what + " should be " + std::to_string(want) + " (got " + std::to_string(got) + ")");
}

// Current value of one observable, by name.
double obs(rulemonkey::RuleMonkeySimulator& sim, const std::string& name) {
  auto names = sim.observable_names();
  auto values = sim.get_observable_values();
  for (size_t i = 0; i < names.size(); ++i)
    if (names[i] == name)
      return values[i];
  check(false, "observable not found: " + name);
  return -1.0;
}

const std::vector<double>& series(const rulemonkey::Result& r, const std::string& name) {
  for (size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == name)
      return r.observable_data[i];
  static const std::vector<double> empty;
  check(false, "observable not found: " + name);
  return empty;
}

long census_count(const std::vector<rulemonkey::SpeciesRow>& rows, const std::string& species) {
  for (const auto& row : rows)
    if (row.species == species)
      return row.count;
  return 0;
}

long census_total(const std::vector<rulemonkey::SpeciesRow>& rows) {
  long total = 0;
  for (const auto& row : rows)
    total += row.count;
  return total;
}

// The model only exercises the component mapping because its seed
// species list `t` before the two `s` components, which is NOT the
// molecule type's declaration order and NOT what BNG2's writeXML emits.
// Regenerating the XML from the .bngl would restore declaration order,
// make the mapping the identity, and quietly cost this file most of what
// it checks — so fail loudly instead.
void test_model_still_lists_components_out_of_order(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator const sim(xml);
  auto rows = sim.initial_species();
  check(rows.size() == 4,
        "the model should declare 4 seed species (got " + std::to_string(rows.size()) + ")");
  if (rows.size() < 2)
    return;
  check(rows[0].name == "P(t~b,s~0,s~1)" && rows[1].name == "P(t~a!1,s~1,s~1).Q(u~y!1)",
        "the seed species must still list `t` before the two `s` components — the "
        "component mapping this model exists to exercise is the identity otherwise.  "
        "seed_build_model.xml is hand-edited on top of BNG2's output for exactly this; "
        "do not regenerate it (got \"" +
            rows[0].name + "\" and \"" + rows[1].name + "\")");
}

// Every copy of every seed species, built from a template that was
// resolved before the first copy existed.
void test_seed_lands_on_every_copy(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(3);

  check_eq(obs(sim, "Pss01"), 3.0, "P monomers with one s~0 and one s~1");
  check_eq(obs(sim, "Pss1"), 4.0, "dimer P's with both s in state 1");
  check_eq(obs(sim, "Ptb"), 3.0, "P monomers with t~b free");
  check_eq(obs(sim, "PtaFree"), 0.0, "P's with t~a FREE — the dimer's t is bonded");
  check_eq(obs(sim, "PtaBound"), 4.0, "P's with t~a bonded");
  check_eq(obs(sim, "Bound"), 4.0, "P(t~a!1).Q(u~y!1) dimers");
  check_eq(obs(sim, "QyBound"), 4.0, "Q's with u~y bonded");
  check_eq(obs(sim, "Qx"), 0.0, "Q(u~x): its seed species has a count of zero");
  check_eq(obs(sim, "Pn"), 7.0, "P molecules");
  check_eq(obs(sim, "Qn"), 4.0, "Q molecules");
  check_eq(obs(sim, "Rn"), 0.0, "R molecules");

  check(sim.get_molecule_count("P") == 7, "the pool should hold 7 P");
  check(sim.get_molecule_count("Q") == 4, "the pool should hold 4 Q");
  check(sim.get_molecule_count("R") == 0, "the pool should hold no R");

  // The structural reading of the same claim: if any copy differed from
  // any other, or from what the species says, the census would carry an
  // extra row rather than two rows of 3 and 4.
  auto rows = sim.enumerate_species();
  check(rows.size() == 2, "the seed pool should hold exactly 2 distinct species (got " +
                              std::to_string(rows.size()) + ")");
  check_eq(static_cast<double>(census_count(rows, "P(s~0,s~1,t~b)")), 3.0,
           "census count of P(s~0,s~1,t~b)");
  check_eq(static_cast<double>(census_count(rows, "P(s~1,s~1,t~a!1).Q(u~y!1)")), 4.0,
           "census count of P(s~1,s~1,t~a!1).Q(u~y!1)");
  sim.destroy_session();
}

// The copy counts drive both the walk and the pool's capacity hint, so
// move all three of them — including the one that was zero, whose seed
// species the walk now skips outright.
void test_counts_scale_with_the_seed_amounts(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_param("nP", 50);
  sim.set_param("nD", 6);
  sim.set_param("nQ", 5);
  sim.initialize(3);

  check_eq(obs(sim, "Pss01"), 50.0, "P monomers after nP -> 50");
  check_eq(obs(sim, "Pss1"), 6.0, "dimer P's after nD -> 6");
  check_eq(obs(sim, "Bound"), 6.0, "dimers after nD -> 6");
  check_eq(obs(sim, "Qx"), 5.0, "free Q(u~x) after nQ 0 -> 5");
  check_eq(obs(sim, "Pn"), 56.0, "P molecules (50 monomers + 6 in dimers)");
  check_eq(obs(sim, "Qn"), 11.0, "Q molecules (6 in dimers + 5 free)");

  auto rows = sim.enumerate_species();
  check(rows.size() == 3,
        "three seed species now have copies (got " + std::to_string(rows.size()) + ")");
  check_eq(static_cast<double>(census_count(rows, "P(s~0,s~1,t~b)")), 50.0,
           "census count of P(s~0,s~1,t~b)");
  check_eq(static_cast<double>(census_count(rows, "P(s~1,s~1,t~a!1).Q(u~y!1)")), 6.0,
           "census count of P(s~1,s~1,t~a!1).Q(u~y!1)");
  check_eq(static_cast<double>(census_count(rows, "Q(u~x)")), 5.0, "census count of Q(u~x)");
  sim.destroy_session();
}

// Past the seed.  No rule creates or destroys a P or a Q, and none of
// them merges two complexes, so every P sits in its own complex for the
// whole run and each is in exactly one of the two `s` shapes — which
// makes Pss01 + Pss1 = 7 an identity, not a statistic.  Likewise each of
// the four t~a slots is either bonded or free, and the only thing it can
// be bonded to is a Q(u~y).
//
// The census is taken twice: once before the run, which populates the
// canonical-label cache under the Debug/ASan self-check, and once after,
// where a birth, a split or a state edit that failed to invalidate a
// cached label would be caught.
void test_run_past_the_seed(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(17);
  auto const before = sim.enumerate_species();
  check(before.size() == 2, "the pre-run census should still see 2 species");

  auto r = sim.simulate(0.0, 20.0, 20);
  const auto& pss01 = series(r, "Pss01");
  const auto& pss1 = series(r, "Pss1");
  const auto& ptafree = series(r, "PtaFree");
  const auto& ptabound = series(r, "PtaBound");
  const auto& bound = series(r, "Bound");
  const auto& qybound = series(r, "QyBound");
  const auto& pn = series(r, "Pn");
  const auto& qn = series(r, "Qn");
  const auto& rn = series(r, "Rn");

  for (size_t i = 0; i < pn.size(); ++i) {
    std::string const at = "[" + std::to_string(i) + "]";
    check_eq(pn[i], 7.0, "P molecules are conserved" + at);
    check_eq(qn[i], 4.0, "Q molecules are conserved" + at);
    check_eq(pss01[i] + pss1[i], 7.0, "every P is in exactly one s shape" + at);
    check_eq(ptafree[i] + ptabound[i], 4.0, "every t~a is bonded or free" + at);
    check_eq(bound[i], ptabound[i], "a bonded t~a is bonded to a Q" + at);
    check_eq(qybound[i], bound[i], "a bonded u~y is bonded to a P" + at);
  }
  // Otherwise the arm proves nothing about mutation: it has to have
  // fired the split, the state flip and the synthesis.
  check(bound.back() < 4.0, "the unbinding rule should have fired (test is otherwise vacuous)");
  check(pss1.back() > 4.0, "the state-flip rule should have fired (test is otherwise vacuous)");
  check(rn.back() > 0.0, "the synthesis rule should have fired (test is otherwise vacuous)");

  // Every complex in the pool, counted two ways.  Seven hold a P; each
  // Q that is not bonded to one is a complex of its own; each
  // synthesized R is a complex of its own.
  auto const after = sim.enumerate_species();
  double const expected = 7.0 + (4.0 - bound.back()) + rn.back();
  check_eq(static_cast<double>(census_total(after)), expected,
           "the post-run census should account for every complex");
  check(after.size() > before.size(),
        "the run should have produced species the seed did not have (got " +
            std::to_string(after.size()) + " vs " + std::to_string(before.size()) + ")");
  sim.destroy_session();
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <seed_build_model.xml>\n", argv[0]);
    return 2;
  }
  const std::string xml = argv[1];

  try {
    test_model_still_lists_components_out_of_order(xml);
    test_seed_lands_on_every_copy(xml);
    test_counts_scale_with_the_seed_amounts(xml);
    test_run_past_the_seed(xml);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: every copy of every seed species is built from its species, and the "
                       "pool stays consistent once the session runs past the seed\n");
  return 0;
}
