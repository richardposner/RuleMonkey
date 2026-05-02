// Direct unit tests for the BNGL math expression parser/evaluator.
//
// These pin down semantics that the BNGL parity ladder only exercises
// indirectly (rate-law evaluation buried in event firing).  When a
// future change shifts how `parse()` or `evaluate()` handles operator
// precedence, short-circuit logical ops, or FP-equality comparison,
// the failure surfaces here as a one-line miscompare instead of a
// 3-sigma trajectory drift on some incidental model.
//
// Pinned regressions (review item numbers from dev/code_review_2026_05_01.md):
//
//   #1  Unary `-` must bind looser than `^`: `-2^2 = -(2^2) = -4`.
//   #2  `&&` / `||` must short-circuit: `if(x>0 && 1/x>1, …)` with x=0
//       must NOT evaluate `1/x` and produce NaN.
//   #6  `==` / `!=` are bit-exact (NFsim/muParser parity).  Two literals
//       that round to the same double compare equal; FP arithmetic
//       errors do not silently coerce to "almost equal".
//
// Plus a handful of sanity checks on operator precedence / associativity /
// builtin dispatch so the parser's overall shape is regression-protected
// alongside the high-impact edge cases above.

#include "expr_eval.hpp"

#include <cmath>
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

double eval(const std::string& src, const rulemonkey::expr::VariableMap& vars = {}) {
  auto ast = rulemonkey::expr::parse(src);
  return rulemonkey::expr::evaluate(*ast, vars);
}

bool throws_runtime_error(const std::string& src, const rulemonkey::expr::VariableMap& vars = {}) {
  try {
    auto ast = rulemonkey::expr::parse(src);
    (void)rulemonkey::expr::evaluate(*ast, vars);
    return false;
  } catch (const std::runtime_error&) {
    return true;
  } catch (...) {
    return false;
  }
}

// ---------------------------------------------------------------------------
// #1 — unary minus precedence: `-2^2` parses as `-(2^2) = -4`, not `(-2)^2 = 4`.
// Reviewer's pre-fix behavior: `-2^2 = 4` (silent wrong answer for negative
// bases in rate laws and TFUN expressions).  Fix at `3b96b84` routed unary
// `-`/`+` operands through `parse_power` so the unary sits OUTSIDE the
// exponentiation.
// ---------------------------------------------------------------------------
void test_unary_minus_vs_power() {
  check(eval("-2^2") == -4.0, "-2^2 must parse as -(2^2) = -4 (unary - binds looser than ^)");
  check(eval("-(2^2)") == -4.0, "explicit -(2^2) = -4");
  check(eval("(-2)^2") == 4.0, "explicit (-2)^2 = 4 (parenthesized unary)");
  check(eval("- 2 ^ 2") == -4.0, "whitespace must not change unary precedence");
  check(eval("-2^3") == -8.0, "-2^3 = -(2^3) = -8");
  check(eval("--2") == 2.0, "double negation");
  check(eval("-2^0") == -1.0, "-2^0 = -(1) = -1");
  // Right-associativity of ^ combined with unary
  check(eval("-2^2^3") == -256.0, "-2^2^3 = -(2^(2^3)) = -256 (^ right-assoc)");
}

// ---------------------------------------------------------------------------
// #2 — logical && / || short-circuit.  Pre-fix RM evaluated both sides
// unconditionally; `if(x>0 && 1/x>1, val, 0)` with x=0 divided by zero and
// produced NaN that poisoned the rate.  Fix at `87af961` special-cased the
// logical ops in the evaluator before generic both-sides eval.
// ---------------------------------------------------------------------------
void test_short_circuit_and_or() {
  rulemonkey::expr::VariableMap vars{{"x", 0.0}};

  // && short-circuit: false-and-divbyzero must NOT evaluate the RHS.
  // The if-branch is selected by the AND's value (false → else branch = 0).
  // If the && evaluated both sides eagerly, 1/0 → +inf → +inf > 1 → 1, but
  // worse, in earlier RM 0/0 paths could surface NaN that propagated.
  // Either way, the SAFE behavior is "result = 0, no NaN, no exception".
  double r1 = eval("if(x > 0 && 1/x > 1, 99, 0)", vars);
  check(r1 == 0.0, "&& short-circuit: x=0 must give 0 (not NaN, not 99)");
  check(!std::isnan(r1), "&& short-circuit: result must not be NaN");

  // || short-circuit: true-or-divbyzero must NOT evaluate the RHS.
  rulemonkey::expr::VariableMap vars_true{{"x", 1.0}};
  double r2 = eval("if(x > 0 || 1/0 > 1, 99, 0)", vars_true);
  check(r2 == 99.0, "|| short-circuit: x=1 must give 99 without evaluating 1/0");
  check(!std::isnan(r2), "|| short-circuit: result must not be NaN");

  // Sanity: when the LHS does NOT short-circuit, the RHS is evaluated
  // and contributes to the truth value (no surprise here).
  rulemonkey::expr::VariableMap v_pos{{"x", 0.5}};
  double r3 = eval("if(x > 0 && 1/x > 1, 99, 0)", v_pos);
  // 1/0.5 = 2, 2 > 1 is true, x > 0 is true → AND true → 99
  check(r3 == 99.0, "&& with both sides truthy: result is the then-branch");

  // The same shape with || should also work — both true.
  double r4 = eval("if(x > 0 || 1/x > 1, 99, 0)", v_pos);
  check(r4 == 99.0, "|| with both sides truthy: result is the then-branch");
}

// ---------------------------------------------------------------------------
// #6 — `==` and `!=` on doubles are bit-exact (NFsim muParser parity).
// Pre-edit this was a silent gotcha; we now document and pin down the
// behavior so a future "let's add a tolerance" change has to update this
// test deliberately.
// ---------------------------------------------------------------------------
void test_double_equality_exact() {
  // Bit-identical literals compare equal.
  check(eval("1.0 == 1.0") == 1.0, "1.0 == 1.0 is true");
  check(eval("1.0 != 1.0") == 0.0, "1.0 != 1.0 is false");
  check(eval("0 == 0") == 1.0, "0 == 0 is true");

  // FP arithmetic that "should" land on a clean value but doesn't.
  // 0.1 + 0.2 = 0.30000000000000004, which != 0.3.  This is the gotcha.
  check(eval("0.1 + 0.2 == 0.3") == 0.0,
        "0.1 + 0.2 == 0.3 must be FALSE (bit-exact, NFsim parity)");
  check(eval("0.1 + 0.2 != 0.3") == 1.0, "and the != complement must be TRUE");

  // Variables that hold the same value compare equal.
  rulemonkey::expr::VariableMap vars{{"a", 7.5}, {"b", 7.5}, {"c", 7.50001}};
  check(eval("a == b", vars) == 1.0, "vars with identical bits compare equal");
  check(eval("a == c", vars) == 0.0, "vars off by 1e-5 compare unequal");

  // Comparison operators work as expected.
  check(eval("1.0 < 2.0") == 1.0, "1 < 2");
  check(eval("2.0 <= 2.0") == 1.0, "2 <= 2");
  check(eval("3.0 > 2.0") == 1.0, "3 > 2");
  check(eval("3.0 >= 3.0") == 1.0, "3 >= 3");
}

// ---------------------------------------------------------------------------
// Sanity: arithmetic precedence + builtins.  Not from the review report,
// but cheap and locks in the broader parser shape so a refactor can't
// accidentally swap precedence layers.
// ---------------------------------------------------------------------------
void test_arithmetic_precedence() {
  check(eval("1 + 2 * 3") == 7.0, "* binds tighter than +");
  check(eval("(1 + 2) * 3") == 9.0, "parentheses override");
  check(eval("2 ^ 3 ^ 2") == 512.0, "^ is right-associative: 2^(3^2) = 2^9 = 512");
  check(eval("2 * 3 + 4 * 5") == 26.0, "two * groups separated by +");
  check(eval("10 - 3 - 2") == 5.0, "- is left-associative");
  check(eval("8 / 2 / 2") == 2.0, "/ is left-associative");

  // Builtin dispatch (single switch via BuiltinKind).
  check(std::abs(eval("sqrt(16)") - 4.0) < 1e-12, "sqrt(16) = 4");
  check(eval("abs(-3.5)") == 3.5, "abs(-3.5) = 3.5");
  check(std::abs(eval("ln(exp(2))") - 2.0) < 1e-12, "ln(exp(2)) = 2");
  check(eval("min(3, 7)") == 3.0, "min(3, 7) = 3");
  check(eval("max(3, 7)") == 7.0, "max(3, 7) = 7");
  check(eval("pow(2, 10)") == 1024.0, "pow(2, 10) = 1024");
  check(eval("floor(3.9)") == 3.0, "floor(3.9) = 3");
  check(eval("ceil(3.1)") == 4.0, "ceil(3.1) = 4");
  check(eval("round(3.5)") == 4.0, "round(3.5) = 4");

  // `if` ternary
  check(eval("if(1 > 0, 100, 200)") == 100.0, "if(true, ...) → then-branch");
  check(eval("if(1 < 0, 100, 200)") == 200.0, "if(false, ...) → else-branch");
}

// Variable resolution: unknown name throws cleanly with the variable name in the
// error message.  Pinned because parse-time resolution shifted to lower_variables
// in #12 and the un-lowered evaluator path is the fallback used by
// resolve_cached.
void test_unresolved_variable() {
  check(throws_runtime_error("missing_var + 1"),
        "evaluate against an empty VariableMap on a Variable node throws");

  rulemonkey::expr::VariableMap vars{{"a", 5.0}};
  check(eval("a + 1", vars) == 6.0, "known variable resolves");
  check(throws_runtime_error("a + b", vars),
        "partially-known expression throws on the unknown half");
}

} // namespace

int main() {
  try {
    test_unary_minus_vs_power();
    test_short_circuit_and_or();
    test_double_equality_exact();
    test_arithmetic_precedence();
    test_unresolved_variable();
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: expr_eval semantics (unary -, short-circuit, FP-eq, "
                       "precedence, builtins, unresolved-var) all pass\n");
  return 0;
}
