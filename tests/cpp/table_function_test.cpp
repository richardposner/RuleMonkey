// Direct unit tests for TableFunction (TFUN) interpolation.
//
// The BNGL parity ladder exercises tfun via models that consume rates
// derived from a table (see ft_tfun.bngl), but the interpolation
// semantics are easier to pin down here in isolation — the parity test
// can pass while the boundary / step semantics drift if the model
// happens to never sample at a breakpoint.
//
// Pinned regressions (review item numbers from dev/code_review_2026_05_01.md):
//
//   #7  Step interpolation is RIGHT-continuous (matches BNG Network3):
//       f(x) = ys[i] for x ∈ [xs[i], xs[i+1]).  At the breakpoint
//       xs[i+1] the value JUMPS to ys[i+1] from the right.  RM's prior
//       comment said "left-continuous", which was wrong terminology
//       (BNG itself mislabels it the same way).  Implementation matches
//       BNG exactly; we test that and lock it down.
//
// Boundary clamp behavior is also pinned (RM matches BNG Network3:
// `if (x <= xs.front()) return ys.front();`, similarly at the upper
// boundary), distinct from NFsim's TFUN cursor which returns 0 before
// the first sample.  Authors who want zero-before-table behavior must
// place an explicit `(xs[0], 0)` leading sample.

#include "table_function.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
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

bool eq(double a, double b, double eps = 1e-12) { return std::abs(a - b) <= eps; }

// ---------------------------------------------------------------------------
// #7 — step interpolation is right-continuous: f(xs[i]) = ys[i],
// f(midpoint of [xs[i], xs[i+1])) = ys[i], and f(xs[i+1]) = ys[i+1]
// (the value JUMPS at the breakpoint from the right).
// ---------------------------------------------------------------------------
void test_step_right_continuous() {
  // Three breakpoints at 1.0 / 2.0 / 3.0 with values 10 / 20 / 30.
  rulemonkey::TableFunction const tfun("step_test", {1.0, 2.0, 3.0}, {10.0, 20.0, 30.0}, "t",
                                       rulemonkey::TfunMethod::Step);

  // At each breakpoint exactly: value is the y-value of the interval STARTING
  // here.  At xs[i+1] (=2.0), f == ys[i+1] (=20).  This is the right-
  // continuous semantic.
  check(eq(tfun.evaluate(1.0), 10.0), "step: f(xs[0]=1.0) = ys[0] = 10");
  check(eq(tfun.evaluate(2.0), 20.0), "step: f(xs[1]=2.0) = ys[1] = 20 (jumps from 10 to 20)");
  check(eq(tfun.evaluate(3.0), 30.0), "step: f(xs[2]=3.0) = ys[2] = 30 (clamp at upper bound)");

  // Midpoints of [xs[i], xs[i+1]) take ys[i] (the LEFT y-value of the bracket).
  check(eq(tfun.evaluate(1.5), 10.0), "step: f(1.5) = ys[0] = 10 (mid of [1,2))");
  check(eq(tfun.evaluate(1.999), 10.0), "step: f(1.999) = ys[0] = 10 (just below xs[1])");
  check(eq(tfun.evaluate(2.5), 20.0), "step: f(2.5) = ys[1] = 20 (mid of [2,3))");
  check(eq(tfun.evaluate(2.999), 20.0), "step: f(2.999) = ys[1] = 20 (just below xs[2])");
}

// ---------------------------------------------------------------------------
// Linear interpolation: f(midpoint) = average of endpoint y values.
// ---------------------------------------------------------------------------
void test_linear_interpolation() {
  rulemonkey::TableFunction const tfun("lin_test", {0.0, 1.0, 2.0}, {0.0, 10.0, 30.0}, "t",
                                       rulemonkey::TfunMethod::Linear);

  // Exact breakpoints
  check(eq(tfun.evaluate(0.0), 0.0), "linear: f(xs[0]) = ys[0]");
  check(eq(tfun.evaluate(1.0), 10.0), "linear: f(xs[1]) = ys[1]");
  check(eq(tfun.evaluate(2.0), 30.0), "linear: f(xs[2]) = ys[2]");

  // Midpoints — linear so mid is mean of bracket values.
  check(eq(tfun.evaluate(0.5), 5.0), "linear: f(0.5) = (0+10)/2 = 5");
  check(eq(tfun.evaluate(1.5), 20.0), "linear: f(1.5) = (10+30)/2 = 20");

  // Quarter and three-quarter points
  check(eq(tfun.evaluate(0.25), 2.5), "linear: f(0.25) = 2.5");
  check(eq(tfun.evaluate(1.75), 25.0), "linear: f(1.75) = 25");
}

// ---------------------------------------------------------------------------
// Boundary clamp: out-of-range x clamps to the nearest endpoint y.
// Matches BNG Network3.  NFsim's TFUN cursor returns 0 before the first
// sample (different convention), but RM and BNG agree on "clamp."
// ---------------------------------------------------------------------------
void test_boundary_clamp() {
  rulemonkey::TableFunction const lin("lin_clamp", {1.0, 2.0, 3.0}, {10.0, 20.0, 30.0}, "t",
                                      rulemonkey::TfunMethod::Linear);
  rulemonkey::TableFunction const step("step_clamp", {1.0, 2.0, 3.0}, {10.0, 20.0, 30.0}, "t",
                                       rulemonkey::TfunMethod::Step);

  // Below xs.front() — clamp to ys.front()
  check(eq(lin.evaluate(0.5), 10.0), "linear: f(below xs.front()) clamps to ys.front()");
  check(eq(lin.evaluate(-1e6), 10.0), "linear: extreme negative clamps to ys.front()");
  check(eq(step.evaluate(0.5), 10.0), "step: f(below xs.front()) clamps to ys.front()");
  check(eq(step.evaluate(-1e6), 10.0), "step: extreme negative clamps to ys.front()");

  // Above xs.back() — clamp to ys.back()
  check(eq(lin.evaluate(3.5), 30.0), "linear: f(above xs.back()) clamps to ys.back()");
  check(eq(lin.evaluate(1e6), 30.0), "linear: extreme positive clamps to ys.back()");
  check(eq(step.evaluate(3.5), 30.0), "step: f(above xs.back()) clamps to ys.back()");
  check(eq(step.evaluate(1e6), 30.0), "step: extreme positive clamps to ys.back()");

  // Exactly at the boundary — covered by the exact-breakpoint cases above
  // (same code path: `x <= xs.front()` and `x >= xs.back()` both fire on
  // equality).  Pin it explicitly so the equality case can't drift to a
  // strict comparison without breaking the test.
  check(eq(lin.evaluate(1.0), 10.0), "linear: f(xs.front()) clamps to ys.front()");
  check(eq(lin.evaluate(3.0), 30.0), "linear: f(xs.back()) clamps to ys.back()");
}

// ---------------------------------------------------------------------------
// Constructor validation: at least 2 points, equal-length x/y, strictly
// monotonic xs.  Pin the contract so a future loosening accidentally
// allows malformed tables to construct.
// ---------------------------------------------------------------------------
void test_constructor_validation() {
  auto throws = [](auto&& fn) {
    try {
      fn();
    } catch (const std::runtime_error&) {
      return true;
    } catch (...) { // NOLINT(bugprone-empty-catch)
    }
    return false;
  };

  check(throws([] { rulemonkey::TableFunction("too_short", {1.0}, {10.0}, "t"); }),
        "ctor: requires at least 2 data points");

  check(throws([] {
          rulemonkey::TableFunction("len_mismatch", {1.0, 2.0}, {10.0}, "t");
        }),
        "ctor: x and y arrays must have equal length");

  check(throws([] {
          rulemonkey::TableFunction("non_monotonic", {2.0, 1.0}, {10.0, 20.0}, "t");
        }),
        "ctor: xs must be strictly increasing (not decreasing)");

  check(throws([] {
          rulemonkey::TableFunction("dup_x", {1.0, 1.0}, {10.0, 20.0}, "t");
        }),
        "ctor: xs must be strictly increasing (no duplicates)");

  // Sanity: the minimal valid table.
  rulemonkey::TableFunction const ok("ok", {0.0, 1.0}, {0.0, 1.0}, "t");
  check(eq(ok.evaluate(0.5), 0.5), "minimal 2-point linear table interpolates");
}

} // namespace

int main() {
  try {
    test_step_right_continuous();
    test_linear_interpolation();
    test_boundary_clamp();
    test_constructor_validation();
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: table_function semantics (step right-continuous, linear interp, "
                       "boundary clamp, ctor validation) all pass\n");
  return 0;
}
