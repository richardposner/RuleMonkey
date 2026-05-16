// Direct unit tests for RuleMonkey's expression layer.
//
// As of issue #6, RM evaluates BNGL rate-law / function / parameter
// expressions with bngsim::ExprTkEvaluator (vendored ExprTk wrapper),
// and keeps only a lightweight dependency scanner — expr::collect_variables
// — in its own code.  This file pins down both:
//
//   1. expr::collect_variables — RM's own code: every genuine identifier
//      is returned, numeric literals (incl. scientific notation) never
//      leak a spurious variable.
//
//   2. ExprTk semantics that the BNGL parity ladder only exercises
//      indirectly.  These were the migration's feared landmines
//      (dev/exprtk_swap_plan_2026_05_16.md §2); pinning them here means a
//      future re-vendor of exprtk.hpp that shifts any of them fails as a
//      one-line miscompare instead of a 3-sigma trajectory drift:
//
//        - unary `-` binds looser than `^`: `-2^2 = -(2^2) = -4`
//        - `^` is right-associative: `2^3^2 = 2^9 = 512`
//        - `==` / `!=` are bit-exact (NFsim/muParser parity), NOT
//          epsilon-based
//        - 3-arg `if(c,a,b)` lazily evaluates only the taken branch
//        - the `ln` alias is natural log

#include "expr_eval.hpp"

#include "bngsim/expression.hpp"

#include <algorithm>
#include <cmath>
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

// Compile + evaluate a standalone expression (no variables).
double eval(const std::string& src) {
  bngsim::ExprTkEvaluator ev;
  int const id = ev.compile(src);
  return ev.evaluate(id);
}

bool collected(const std::string& expr, const std::string& name) {
  auto vars = rulemonkey::expr::collect_variables(expr);
  return std::find(vars.begin(), vars.end(), name) != vars.end();
}

// ---------------------------------------------------------------------------
// expr::collect_variables — RM's dependency scanner.
// ---------------------------------------------------------------------------
void test_collect_variables() {
  // Genuine identifiers are all returned.
  check(collected("k_break*(1+(0.30*A_bound))", "k_break"), "collect: k_break");
  check(collected("k_break*(1+(0.30*A_bound))", "A_bound"), "collect: A_bound");

  // A local-function argument inside a call survives the scan (harmless —
  // it simply matches no parameter/observable/function downstream).
  check(collected("k_base*exp((-J)*NbrUp(z))", "NbrUp"), "collect: function ref NbrUp");
  check(collected("k_base*exp((-J)*NbrUp(z))", "z"), "collect: local arg z");
  check(collected("k_base*exp((-J)*NbrUp(z))", "J"), "collect: J");

  // Numeric literals must NOT leak identifiers — in particular the
  // exponent marker of a scientific-notation literal.
  check(!collected("1e-5 + foo", "e"), "collect: `1e-5` must not yield variable `e`");
  check(collected("1e-5 + foo", "foo"), "collect: foo alongside `1e-5`");
  check(!collected("3.0 * 2", "3"), "collect: digits are not identifiers");

  // De-duplication: a name referenced twice appears once.
  auto vars = rulemonkey::expr::collect_variables("a + a*b + a");
  check(std::count(vars.begin(), vars.end(), std::string("a")) == 1, "collect: dedup `a`");

  // Comparison / conditional expression — every operand identifier shows up.
  check(collected("if((A<threshold),k_on,k_off)", "A"), "collect: A in if-cond");
  check(collected("if((A<threshold),k_on,k_off)", "threshold"), "collect: threshold");
  check(collected("if((A<threshold),k_on,k_off)", "k_on"), "collect: k_on");
  check(collected("if((A<threshold),k_on,k_off)", "k_off"), "collect: k_off");
}

// ---------------------------------------------------------------------------
// Unary minus vs. `^`, and `^` associativity.
// ---------------------------------------------------------------------------
void test_unary_minus_and_power() {
  check(eval("-2^2") == -4.0, "-2^2 must parse as -(2^2) = -4 (unary - looser than ^)");
  check(eval("-(2^2)") == -4.0, "explicit -(2^2) = -4");
  check(eval("(-2)^2") == 4.0, "explicit (-2)^2 = 4");
  check(eval("-2^3") == -8.0, "-2^3 = -(2^3) = -8");
  check(eval("2^3^2") == 512.0, "^ is right-associative: 2^(3^2) = 2^9 = 512");
}

// ---------------------------------------------------------------------------
// `==` / `!=` are bit-exact, not epsilon-based.
// ---------------------------------------------------------------------------
void test_double_equality_exact() {
  check(eval("1.0 == 1.0") == 1.0, "1.0 == 1.0 is true");
  check(eval("1.0 != 1.0") == 0.0, "1.0 != 1.0 is false");

  // The classic FP gotcha: 0.1 + 0.2 = 0.30000000000000004 != 0.3.
  check(eval("0.1 + 0.2 == 0.3") == 0.0, "0.1 + 0.2 == 0.3 must be FALSE (bit-exact)");
  check(eval("0.1 + 0.2 != 0.3") == 1.0, "and the != complement must be TRUE");

  // A sub-epsilon-but-nonzero difference still compares unequal.
  check(eval("1e-11 == 0") == 0.0, "1e-11 == 0 must be FALSE (no epsilon window)");

  check(eval("1.0 < 2.0") == 1.0, "1 < 2");
  check(eval("2.0 <= 2.0") == 1.0, "2 <= 2");
  check(eval("3.0 >= 3.0") == 1.0, "3 >= 3");
}

// ---------------------------------------------------------------------------
// 3-arg `if` lazily evaluates only the taken branch — the untaken branch
// must not run (e.g. a division by zero in it must not surface).
// ---------------------------------------------------------------------------
void test_if_lazy_ternary() {
  double const r1 = eval("if(0, 1/0, 5)");
  check(r1 == 5.0, "if(false, …) selects the else-branch");
  check(!std::isnan(r1), "untaken then-branch (1/0) must not poison the result");

  double const r2 = eval("if(1, 5, 1/0)");
  check(r2 == 5.0, "if(true, …) selects the then-branch");
  check(!std::isnan(r2), "untaken else-branch (1/0) must not poison the result");
}

// ---------------------------------------------------------------------------
// Builtins + the `ln` alias (natural log), and variable binding.
// ---------------------------------------------------------------------------
void test_builtins_and_variables() {
  check(std::abs(eval("sqrt(16)") - 4.0) < 1e-12, "sqrt(16) = 4");
  check(eval("max(3, 7)") == 7.0, "max(3, 7) = 7");
  check(eval("pow(2, 10)") == 1024.0, "pow(2, 10) = 1024");
  check(std::abs(eval("ln(exp(2))") - 2.0) < 1e-12, "ln is natural log: ln(exp(2)) = 2");

  // Variables bound by address — the contract RM's engine relies on:
  // mutate the bound storage, re-evaluate, see the new value.
  bngsim::ExprTkEvaluator ev;
  double a = 3.0;
  double b = 4.0;
  ev.define_variable("a", &a);
  ev.define_variable("b", &b);
  int const id = ev.compile("a*a + b*b");
  check(ev.evaluate(id) == 25.0, "a*a + b*b with a=3,b=4 = 25");
  a = 6.0;
  b = 8.0;
  check(ev.evaluate(id) == 100.0, "re-evaluate after mutating bound storage = 100");
}

} // namespace

int main() {
  try {
    test_collect_variables();
    test_unary_minus_and_power();
    test_double_equality_exact();
    test_if_lazy_ternary();
    test_builtins_and_variables();
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: expr layer (collect_variables, unary -, ^ assoc, "
                       "bit-exact ==, lazy if, builtins, var binding) all pass\n");
  return 0;
}
