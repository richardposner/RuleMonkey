// Arrhenius rules the load-time expansion cannot price exactly are refused
// (issue #109).
//
// RuleMonkey expands each Arrhenius rule at load into conventional rules with
// fixed rates (energy_expand.{hpp,cpp}).  The expansion turns each energy
// pattern holding the rule's bond into at most one "bound to anything"
// condition on one of the two molecules, and builds every rule from the two
// bond sites alone.  Outside the shapes that survive that, the rule used to
// run at the wrong rate with no warning: the issue measured each fixture
// below 48-52 standard errors off BioNetGen's network, or, for the rule's own
// reactant context, firing on complexes BioNetGen never lets react.  Each is
// now a Tier-0 Error naming the rule and the pattern or context that put it
// outside.
//
// Fixtures (tests/cpp/arrhenius_context/, from the issue's reproducers):
//   ep_state                 a state in an energy pattern           (case 1)
//   ep_free_site             a site an energy pattern requires free (case 2)
//   ep_wildcard_only         wildcard bonds only, L(r!+,r!+)        (case 3)
//   ep_wildcard_center_name  a `!+` site named like the reaction    (case 4)
//                            site
//   ep_repeated_site         a context bond on a second site named  (case 5)
//                            like the reaction site
//   ep_partner_ambiguous     a context site with two possible       (case 6)
//                            partners
//   ep_two_bonds_out         context two bonds out                  (case 7)
//   ep_two_context_bonds     two context bonds                      (case 8)
//   rule_reactant_context    a state in the rule's own reactant     (case 9)
//                            pattern
//   ep_partner_from_rule     case 6, the second partner made only by a rule
//   ep_repeated_context_site a context site its molecule type declares twice
//   ep_wildcard_multiplicity a `!?` that can land two ways, L(r!1,r!?).R(l!1)
//                            on L(r,r,r)
//   rule_context_shapes      the other rule-context checks, one rule each,
//                            and a control rule that must expand
// Controls, passed after the fixture directory: models in the shapes the
// expansion prices exactly (the plain pair, and the pair plus one bond to a
// third molecule type whose site binds nothing else), which must load with no
// Arrhenius diagnostic at all.

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

bool contains(const std::string& haystack, const std::string& needle) {
  return haystack.find(needle) != std::string::npos;
}

// The Arrhenius diagnostics a model loads with.
std::vector<rulemonkey::UnsupportedFeature> arrhenius_features(const std::string& xml) {
  std::vector<rulemonkey::UnsupportedFeature> out;
  rulemonkey::RuleMonkeySimulator const sim(xml);
  for (const auto& f : sim.unsupported_features())
    if (f.element == "RateLaw@type=Arrhenius")
      out.push_back(f);
  return out;
}

struct Refusal {
  std::string rule;   // BNGL rule name
  std::string reason; // the phrase that must appear in its diagnostic
};

// `xml` must refuse exactly the rules in `want`, each with an Error whose text
// names the rule and carries its reason.
void expect_refusals(const std::string& label, const std::string& xml,
                     const std::vector<Refusal>& want) {
  std::vector<rulemonkey::UnsupportedFeature> got;
  try {
    got = arrhenius_features(xml);
  } catch (const std::exception& e) {
    check(false, label + ": load threw: " + e.what());
    return;
  }
  check(got.size() == want.size(), label + ": expected " + std::to_string(want.size()) +
                                       " Arrhenius refusal(s), got " + std::to_string(got.size()));
  for (const auto& w : want) {
    bool found = false;
    for (const auto& f : got) {
      if (!contains(f.feature, "rule '" + w.rule + "'"))
        continue;
      found = true;
      check(f.severity == rulemonkey::Severity::Error,
            label + ": rule " + w.rule + " must be refused with an Error, not a Warn");
      check(contains(f.feature, w.reason), label + ": rule " + w.rule +
                                               " refused without the reason \"" + w.reason +
                                               "\": " + f.feature);
    }
    check(found, label + ": rule " + w.rule + " was not refused");
  }
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <arrhenius_context dir> [control.xml ...]\n", argv[0]);
    return 2;
  }
  std::string const dir = std::string(argv[1]) + "/";

  expect_refusals("case 1", dir + "ep_state.xml",
                  {{"Unbind", "energy pattern A(b!1).B(a!1,s~1) requires B.s in state 1"}});
  expect_refusals("case 2", dir + "ep_free_site.xml",
                  {{"Unbind", "energy pattern A(b!1,c).B(a!1) requires A.c free"}});
  expect_refusals("case 3", dir + "ep_wildcard_only.xml",
                  {{"Unbind", "energy pattern L(r!+,r!+) constrains L.r without naming the "
                              "L.r-R.l bond"}});
  expect_refusals("case 4", dir + "ep_wildcard_center_name.xml",
                  {{"Unbind", "energy pattern L(r!1,r!+).R(l!1) names a second L.r besides the "
                              "reaction site"}});
  expect_refusals("case 5", dir + "ep_repeated_site.xml",
                  {{"Unbind", "energy pattern L(r!1,r!2).R(l!1).S(l!2) names a second L.r "
                              "besides the reaction site"}});
  expect_refusals("case 6", dir + "ep_partner_ambiguous.xml",
                  {{"Unbind", "energy pattern A(b!1,c!2).B(a!1).C(a!2) requires A.c bound to "
                              "C.a, but A.c can also bind D.a"}});
  expect_refusals("case 7", dir + "ep_two_bonds_out.xml",
                  {{"Unbind", "energy pattern A(b!1).B(a!1,c!2).C(b!2,d!3).D(c!3) names more "
                              "than one bond besides the A.b-B.a bond"}});
  expect_refusals("case 8", dir + "ep_two_context_bonds.xml",
                  {{"Unbind", "energy pattern A(b!1,c!2,d!3).B(a!1).C(a!2).D(a!3) names more "
                              "than one bond besides the A.b-B.a bond"}});
  expect_refusals("case 9", dir + "rule_reactant_context.xml",
                  {{"Unbind", "its reactant pattern requires A.s in state 1"}});
  expect_refusals("partner from a rule", dir + "ep_partner_from_rule.xml",
                  {{"Unbind", "energy pattern A(b!1,c!2).B(a!1).C(a!2) requires A.c bound to "
                              "C.a, but A.c can also bind D.a"}});
  expect_refusals("`!?` multiplicity", dir + "ep_wildcard_multiplicity.xml",
                  {{"Unbind", "energy pattern L(r!1,r!?).R(l!1) can place its L.r (!?) 2 ways"}});
  expect_refusals("repeated context site", dir + "ep_repeated_context_site.xml",
                  {{"Unbind", "energy pattern A(x!1,c!2).B(y!1).C(a!2) requires A.c bound, but "
                              "A declares more than one c"}});
  // Ok, the control, must expand: expect_refusals counts every refusal.
  expect_refusals("rule context", dir + "rule_context_shapes.xml",
                  {{"Extra", "its reactant pattern names a molecule besides the two it unbinds"},
                   {"Free", "its reactant pattern constrains the bond of A2.c"},
                   {"Bound", "its reactant pattern constrains the bond of A4.c"},
                   {"Dot", "it breaks the bond without writing its products apart"}});

  for (int i = 2; i < argc; ++i)
    expect_refusals(std::string("control ") + argv[i], argv[i], {});

  if (g_failures != 0) {
    std::fprintf(stderr, "arrhenius_context_test: %d check(s) failed\n", g_failures);
    return 1;
  }
  std::printf("arrhenius_context_test: OK\n");
  return 0;
}
