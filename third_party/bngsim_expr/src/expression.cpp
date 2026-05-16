// bngsim/src/expression.cpp — ExprTk-based expression evaluator
//
// ExprTk replaces muParser and provides BNG-compatible aliases.
// Supports underscore-prefixed constants, built-in functions, and time().
//
// Note: BNG2.pl's muParser exposes simulation time as the zero-arg function
// time() and leaves t free as an ordinary identifier. We mirror that here so
// that BNGL models defining a parameter / observable / species named `t`
// (a common counter pattern, e.g. `Molecules t counter()`) load successfully.

// ExprTk is a large header — compile once here.
// Disable some ExprTk features we don't need to speed up compilation.
#define exprtk_disable_string_capabilities
#define exprtk_disable_rtl_io_file
#define exprtk_disable_rtl_vecops
// BNG is case-sensitive for parameter names (e.g., k3 ≠ K3).
// ExprTk defaults to case-insensitive, which silently merges k3/K3
// into the same variable and produces incorrect trajectories.
#define exprtk_disable_caseinsensitivity
#include "exprtk.hpp"

#include "bngsim/expression.hpp"

#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace bngsim {

// ─── Mratio: confluent hypergeometric ratio M(a+1,b+1,z)/M(a,b,z) ───────────
//
// Ported from BNG's Util::Mratio (W. S. Hlavacek, 2018).
// Uses modified Lentz continued fraction method.
static double mratio_impl(double a, double b, double z) {
    // M(a+1,b+1,z)/M(a,b,z) via continued fraction (Gauss).
    // The continued fraction for the ratio of contiguous 1F1 functions:
    //   M(a+1,b+1,z)/M(a,b,z) = 1/(1 - az/b(b+1) / (1 - (a+1)z/(b+1)(b+2) / ...))
    //
    // Evaluated using the modified Lentz algorithm.

    const double tiny = 1.0e-30;
    const double eps = 1.0e-15;
    const int max_iter = 10000;

    // Use the recurrence: M(a+1,b+1,z)/M(a,b,z)
    // via the continued fraction representation.
    // Reference: Abramowitz & Stegun 13.5; DLMF 13.2.
    //
    // The CF coefficients for the ratio r = M(a+1,b+1,z)/M(a,b,z):
    //   r = b/(b - z + ...) with specific CF terms.
    //
    // Using the simpler series ratio approach from BNG:
    // Compute ratio via Lentz method on the CF:
    //   f = b/(b-z + CF)
    // where CF has terms a_n/b_n.

    // Direct implementation matching BNG's Util::Mratio:
    // Use the Gauss continued fraction for the ratio of Kummer functions.
    //
    // The ratio M(a+1,b+1,z)/M(a,b,z) can be expressed as:
    //   b / (b - z * M(a+1,b+2,z)/M(a+1,b+1,z))
    //
    // This leads to the CF:
    //   r = b / (b - z * (a+1)/(b+1) / (1 + z*(b-a)/(b+1)(b+2) / (1 + ...)))
    //
    // We use the standard Lentz algorithm.

    // Actually, let's use the direct series computation which is more robust
    // for the parameter ranges encountered in BNG models.

    // Compute M(a,b,z) and M(a+1,b+1,z) by series and take ratio.
    // M(a,b,z) = sum_{n=0}^{inf} (a)_n / (b)_n * z^n / n!

    double term0 = 1.0; // running term for M(a,b,z)
    double sum0 = 1.0;  // partial sum for M(a,b,z)
    double term1 = 1.0; // running term for M(a+1,b+1,z)
    double sum1 = 1.0;  // partial sum for M(a+1,b+1,z)

    for (int n = 1; n <= max_iter; ++n) {
        double nn = static_cast<double>(n);
        term0 *= (a + nn - 1.0) / (b + nn - 1.0) * z / nn;
        sum0 += term0;
        term1 *= (a + nn) / (b + nn) * z / nn;
        sum1 += term1;

        // Check convergence of both series
        if (std::abs(term0) < eps * std::abs(sum0) && std::abs(term1) < eps * std::abs(sum1)) {
            return sum1 / sum0;
        }
    }

    // If we didn't converge, return best estimate
    return sum1 / sum0;
}

// ─── ExprTk custom function adapters ─────────────────────────────────────────
//
// ExprTk requires inheriting from ifunction for custom functions.

// 3-arg: mratio(a, b, z)
template <typename T> struct MratioFunction : public exprtk::ifunction<T> {
    MratioFunction() : exprtk::ifunction<T>(3) {
        exprtk::ifunction<T>::allow_zero_parameters() = false;
    }
    T operator()(const T &a, const T &b, const T &z) override {
        return static_cast<T>(
            mratio_impl(static_cast<double>(a), static_cast<double>(b), static_cast<double>(z)));
    }
};

// 1-arg aliases for backward compat
template <typename T> struct LnFunction : public exprtk::ifunction<T> {
    LnFunction() : exprtk::ifunction<T>(1) {}
    T operator()(const T &x) override { return std::log(x); }
};

template <typename T> struct RintFunction : public exprtk::ifunction<T> {
    RintFunction() : exprtk::ifunction<T>(1) {}
    T operator()(const T &x) override { return std::round(x); }
};

template <typename T> struct SignFunction : public exprtk::ifunction<T> {
    SignFunction() : exprtk::ifunction<T>(1) {}
    T operator()(const T &x) override { return (x > 0.0) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0); }
};

// 0-arg: time() — reads from a bound double*
template <typename T> struct TimeFunction : public exprtk::ifunction<T> {
    double *time_ptr = nullptr;
    TimeFunction() : exprtk::ifunction<T>(0) {
        exprtk::ifunction<T>::allow_zero_parameters() = true;
    }
    T operator()() override { return time_ptr ? static_cast<T>(*time_ptr) : T(0); }
};

// Adapter for std::function-based custom functions (0-3 args)
template <typename T> struct StdFunc0Adapter : public exprtk::ifunction<T> {
    std::function<double()> fn;
    StdFunc0Adapter(std::function<double()> f) : exprtk::ifunction<T>(0), fn(std::move(f)) {
        exprtk::ifunction<T>::allow_zero_parameters() = true;
    }
    T operator()() override { return static_cast<T>(fn()); }
};

template <typename T> struct StdFunc1Adapter : public exprtk::ifunction<T> {
    std::function<double(double)> fn;
    StdFunc1Adapter(std::function<double(double)> f) : exprtk::ifunction<T>(1), fn(std::move(f)) {}
    T operator()(const T &x) override { return static_cast<T>(fn(x)); }
};

template <typename T> struct StdFunc2Adapter : public exprtk::ifunction<T> {
    std::function<double(double, double)> fn;
    StdFunc2Adapter(std::function<double(double, double)> f)
        : exprtk::ifunction<T>(2), fn(std::move(f)) {}
    T operator()(const T &x, const T &y) override { return static_cast<T>(fn(x, y)); }
};

template <typename T> struct StdFunc3Adapter : public exprtk::ifunction<T> {
    std::function<double(double, double, double)> fn;
    StdFunc3Adapter(std::function<double(double, double, double)> f)
        : exprtk::ifunction<T>(3), fn(std::move(f)) {}
    T operator()(const T &x, const T &y, const T &z) override {
        return static_cast<T>(fn(x, y, z));
    }
};

// ─── ExprTk implementation (pimpl) ───────────────────────────────────────────

// ─── Name remapping ──────────────────────────────────────────────────────────
//
// Two transparent transformations bridge BNG identifier conventions onto
// ExprTk's identifier rules:
//
//   1. Underscore-prefixed names: ExprTk rejects identifiers starting with
//      '_'. BNG built-ins use _pi, _e, _NA, etc. We map "_X" → "u_X" on
//      registration and rewrite the same tokens in compiled expressions.
//      This mapping is unconditional (no built-in ExprTk symbol starts with
//      '_'), so the rewrite is safe to apply token-by-token in expressions.
//
//   2. ExprTk reserved-word collisions: ExprTk rejects add_variable() for
//      names matching its reserved-word/reserved-symbol lists (e.g., a BNG
//      parameter literally named `const`, `true`, `false`). We register the
//      user's variable under "r_<name>" and rewrite references in compiled
//      expressions. Unlike (1) this rewrite is *conditional* — built-in
//      functions like `sin`, `if`, `time` must remain literal tokens — so
//      we only rewrite identifiers that were actually mangled at
//      registration. The mangling map is per-evaluator and lives on Impl.

static bool is_exprtk_reserved(const std::string &name) {
    // Names that, if registered as a user variable, would collide with
    // a name already taken by the symbol table. Two sources:
    //
    //   1. ExprTk's reserved_words[] + reserved_symbols[] (see
    //      third_party/exprtk/exprtk.hpp). Both lists reject add_variable().
    //
    //   2. bngsim-specific function aliases registered in
    //      ExprTkEvaluator::Impl::init_builtins() (`ln`, `rint`, `sign`,
    //      `mratio`, `time`). ExprTk would reject add_variable() on any
    //      of these once the function is registered, but the names are
    //      not on ExprTk's reserved_symbols[] list, so we have to track
    //      them here ourselves. BNG2.pl's parser already rejects
    //      `ln`/`rint`/`mratio`/`time` as parameter names upstream, so
    //      models reaching bngsim via `generate_network` can only hit
    //      the `sign` collision in practice — but we mangle all five
    //      for symmetry and to handle hand-crafted .net inputs.
    //
    // Comparison is exact (case-sensitive) because we build with
    // exprtk_disable_caseinsensitivity, so e.g. `Const` is not reserved.
    static const std::unordered_set<std::string> reserved = {
        // Control-flow / language keywords (ExprTk reserved_words[])
        "assert", "break", "case", "continue", "const", "default", "false", "for", "if", "else",
        "ilike", "in", "like", "and", "nand", "nor", "not", "null", "or", "repeat", "return", "shl",
        "shr", "swap", "switch", "true", "until", "var", "while", "xnor", "xor",
        // Built-in math / utility functions (ExprTk reserved_symbols[])
        "abs", "acos", "acosh", "asin", "asinh", "atan", "atanh", "atan2", "avg", "ceil", "clamp",
        "cos", "cosh", "cot", "csc", "deg2grad", "deg2rad", "equal", "erf", "erfc", "exp", "expm1",
        "floor", "frac", "grad2deg", "hypot", "iclamp", "inrange", "log", "log10", "log2", "logn",
        "log1p", "mand", "max", "min", "mod", "mor", "mul", "ncdf", "not_equal", "pow", "rad2deg",
        "root", "round", "roundn", "sec", "sgn", "sin", "sinc", "sinh", "sqrt", "sum", "tan",
        "tanh", "trunc",
        // bngsim-registered function aliases (init_builtins)
        "ln", "rint", "sign", "mratio", "time"};
    return reserved.find(name) != reserved.end();
}

// Compute the symbol-table key for `name` at registration time.
// Combines the unconditional underscore remap with reserved-word mangling.
static std::string compute_registration_name(const std::string &name) {
    if (!name.empty() && name[0] == '_') {
        return "u_" + name.substr(1);
    }
    if (is_exprtk_reserved(name)) {
        return "r_" + name;
    }
    return name;
}

struct ExprTkEvaluator::Impl {
    using SymbolTable = exprtk::symbol_table<double>;
    using Expression = exprtk::expression<double>;
    using Parser = exprtk::parser<double>;

    SymbolTable symbol_table;
    // Parser is shared across clones to avoid re-constructing the ~100KB
    // template object. The parser is stateless between compile() calls.
    std::shared_ptr<Parser> parser;
    std::vector<Expression> expressions;

    // Cached preprocessed expression strings for efficient clone().
    // Indexed in parallel with `expressions`.
    std::vector<std::string> preprocessed_strings;

    // Owned custom function objects (must outlive the symbol table)
    MratioFunction<double> mratio_func;
    LnFunction<double> ln_func;
    RintFunction<double> rint_func;
    SignFunction<double> sign_func;
    TimeFunction<double> time_func;

    // User-registered custom functions (owned, heap-allocated)
    std::vector<std::unique_ptr<exprtk::ifunction<double>>> user_functions;

    // Names that were mangled at registration to avoid ExprTk reserved-word
    // collisions (key: original BNG name, value: ExprTk symbol-table key).
    // Underscore-prefixed names (e.g., _pi) are NOT recorded here — they
    // remap unconditionally via compute_registration_name() in both
    // directions, so no per-evaluator state is needed.
    std::unordered_map<std::string, std::string> mangled_user_names;

    // BNG-source names registered as scalar variables/constants on this
    // evaluator (i.e., everything that went through define_variable /
    // add_remapped_constant). Used by strip_empty_parens() to rewrite
    // `obs()` → `obs` for observable references that the BNG parser
    // accepts (Expression.pm:870-927 — Observable as zero-arg call) but
    // ExprTk's grammar would reject as `obs * ()`.
    std::unordered_set<std::string> scalar_variable_names;

    // Look up the symbol-table key for `name` when rewriting an expression.
    // Mirrors compute_registration_name() but only mangles reserved words
    // that were actually registered on this evaluator, so built-in tokens
    // (sin, if, time, ...) pass through unchanged.
    std::string remap_token(const std::string &name) const {
        if (!name.empty() && name[0] == '_') {
            return "u_" + name.substr(1);
        }
        auto it = mangled_user_names.find(name);
        if (it != mangled_user_names.end()) {
            return it->second;
        }
        return name;
    }

    // Rewrite `name()` → `name` for any identifier registered as a scalar
    // variable on this evaluator. BNGL's grammar (per BNG2.pl's
    // bng2/Perl2/Expression.pm:870-927) accepts an Observable as a
    // zero-arg call (`obs()`), and BNG2.pl preserves that syntax verbatim
    // when emitting .net rate laws / parameter expressions / event
    // expressions. ExprTk's grammar would parse `obs()` as the implicit
    // multiplication `obs * ()` and reject the empty parens with ERR248.
    // We close the gap by stripping the trailing `()` for any name we
    // know is a scalar — leaving function names (built-ins like `sin`,
    // `time`, user-defined Func0/1/2/3) untouched, since those go
    // through add_function and are not in scalar_variable_names.
    std::string strip_empty_parens(const std::string &expr) const {
        if (scalar_variable_names.empty())
            return expr;
        std::string result;
        result.reserve(expr.size());
        size_t i = 0;
        while (i < expr.size()) {
            const bool at_boundary =
                (i == 0) ||
                (!std::isalnum(static_cast<unsigned char>(expr[i - 1])) && expr[i - 1] != '_');
            const bool ident_start =
                (std::isalpha(static_cast<unsigned char>(expr[i])) || expr[i] == '_');
            if (at_boundary && ident_start) {
                size_t start = i;
                i++;
                while (i < expr.size() &&
                       (std::isalnum(static_cast<unsigned char>(expr[i])) || expr[i] == '_')) {
                    i++;
                }
                std::string ident = expr.substr(start, i - start);
                if (i + 1 < expr.size() && expr[i] == '(' && expr[i + 1] == ')' &&
                    scalar_variable_names.count(ident)) {
                    result += ident;
                    i += 2;
                } else {
                    result += ident;
                }
                continue;
            }
            result += expr[i];
            i++;
        }
        return result;
    }

    // Rewrite all identifier tokens in `expr` through remap_token().
    std::string remap_expression(const std::string &expr) const {
        std::string result;
        result.reserve(expr.size() + 16);
        size_t i = 0;
        while (i < expr.size()) {
            const bool at_boundary =
                (i == 0) ||
                (!std::isalnum(static_cast<unsigned char>(expr[i - 1])) && expr[i - 1] != '_');
            const bool ident_start =
                (std::isalpha(static_cast<unsigned char>(expr[i])) || expr[i] == '_');
            if (at_boundary && ident_start) {
                size_t start = i;
                i++;
                while (i < expr.size() &&
                       (std::isalnum(static_cast<unsigned char>(expr[i])) || expr[i] == '_')) {
                    i++;
                }
                result += remap_token(expr.substr(start, i - start));
                continue;
            }
            result += expr[i];
            i++;
        }
        return result;
    }

    void add_remapped_constant(const std::string &name, double value) {
        std::string mapped = compute_registration_name(name);
        if (mapped != name && name[0] != '_') {
            mangled_user_names[name] = mapped;
        }
        symbol_table.add_constant(mapped, value);
        scalar_variable_names.insert(name);
    }

    void init_builtins() {
        // Register built-in constants
        // Note: ExprTk rejects '_' prefix, so we remap to "u_" internally
        add_remapped_constant("_pi", 3.14159265358979323846);
        add_remapped_constant("_e", 2.71828182845904523536);
        add_remapped_constant("_NA", 6.02214076e23);
        add_remapped_constant("_kB", 1.380649e-23);
        add_remapped_constant("_R", 8.314462618153241);
        add_remapped_constant("_h", 6.62607015e-34);
        add_remapped_constant("_F", 96485.33212331002);

        // Register backward-compatible aliases
        symbol_table.add_function("ln", ln_func);
        symbol_table.add_function("rint", rint_func);
        symbol_table.add_function("sign", sign_func);

        // Register built-in functions
        // Note: ExprTk has a built-in `if` keyword (grammar-level) that handles
        // if(cond, true_val, false_val) natively. No custom function needed.
        symbol_table.add_function("mratio", mratio_func);

        // time() — pointer set later via set_time_ptr().
        // We do NOT register `t` here so that `t` is free for use as a model
        // identifier (matches BNG2.pl convention).
        symbol_table.add_function("time", time_func);
    }

    Impl() : parser(std::make_shared<Parser>()) {
        // Increase max stack depth for deeply nested if() expressions.
        // ExprTk default is 400 (~200 nested if()), muParser handled 2000.
        parser->settings().set_max_stack_depth(4096);
        init_builtins();
    }

    // Constructor that shares an existing parser (for clone_empty)
    explicit Impl(std::shared_ptr<Parser> shared_parser) : parser(std::move(shared_parser)) {
        init_builtins();
    }

    void set_time_ptr(double *ptr) { time_func.time_ptr = ptr; }
};

ExprTkEvaluator::ExprTkEvaluator() : impl_(std::make_unique<Impl>()) {}

ExprTkEvaluator::~ExprTkEvaluator() = default;

void ExprTkEvaluator::define_variable(const std::string &name, double *addr) {
    std::string mapped = compute_registration_name(name);
    if (!impl_->symbol_table.add_variable(mapped, *addr)) {
        const bool reserved = is_exprtk_reserved(name);
        std::string detail;
        if (reserved) {
            // The mangled form already exists — most often because another
            // BNG name collides with the mangled key (e.g., user has both
            // `const` and a separate `r_const`).
            detail = "name '" + name + "' is an ExprTk reserved word; bngsim mangles it to '" +
                     mapped +
                     "', and that mangled key is already registered. Rename one of the "
                     "conflicting parameters.";
        } else {
            detail = "name '" + name + "' (mapped: '" + mapped +
                     "') is already registered. Check the .net file for duplicate "
                     "parameter / observable / species names (case-sensitive).";
        }
        throw std::runtime_error("ExprTk: failed to register variable '" + name + "'. " + detail);
    }
    if (mapped != name && (name.empty() || name[0] != '_')) {
        impl_->mangled_user_names[name] = mapped;
    }
    impl_->scalar_variable_names.insert(name);
}

void ExprTkEvaluator::define_constant(const std::string &name, double value) {
    impl_->add_remapped_constant(name, value);
}

void ExprTkEvaluator::define_function(const std::string &name, Func0 fn) {
    auto adapter = std::make_unique<StdFunc0Adapter<double>>(std::move(fn));
    impl_->symbol_table.add_function(name, *adapter);
    impl_->user_functions.push_back(std::move(adapter));
}

void ExprTkEvaluator::define_function(const std::string &name, Func1 fn) {
    auto adapter = std::make_unique<StdFunc1Adapter<double>>(std::move(fn));
    impl_->symbol_table.add_function(name, *adapter);
    impl_->user_functions.push_back(std::move(adapter));
}

void ExprTkEvaluator::define_function(const std::string &name, Func2 fn) {
    auto adapter = std::make_unique<StdFunc2Adapter<double>>(std::move(fn));
    impl_->symbol_table.add_function(name, *adapter);
    impl_->user_functions.push_back(std::move(adapter));
}

void ExprTkEvaluator::define_function(const std::string &name, Func3 fn) {
    auto adapter = std::make_unique<StdFunc3Adapter<double>>(std::move(fn));
    impl_->symbol_table.add_function(name, *adapter);
    impl_->user_functions.push_back(std::move(adapter));
}

// ─── C-style logical operator replacement ────────────────────────────────────
// BNG2.pl emits C-style && and || in if() conditions, but ExprTk uses
// 'and' and 'or' keywords. Replace before compilation.
static std::string replace_logical_operators(const std::string &expr) {
    std::string result;
    result.reserve(expr.size() + 16);
    for (size_t i = 0; i < expr.size(); ++i) {
        if (i + 1 < expr.size()) {
            if (expr[i] == '&' && expr[i + 1] == '&') {
                result += " and ";
                i++; // skip second '&'
                continue;
            }
            if (expr[i] == '|' && expr[i + 1] == '|') {
                result += " or ";
                i++; // skip second '|'
                continue;
            }
        }
        result += expr[i];
    }
    return result;
}

int ExprTkEvaluator::compile(const std::string &expr) {
    // Replace C-style logical operators before any other processing
    std::string preprocessed = replace_logical_operators(expr);

    // Strip `obs()` → `obs` for any name registered as a scalar variable.
    // BNGL accepts Observable references as zero-arg calls; ExprTk does
    // not. Run before remap_expression so we match against BNG-source
    // names, not their post-mangling forms.
    std::string stripped = impl_->strip_empty_parens(preprocessed);

    // Remap identifiers before ExprTk compilation:
    //   - Unconditional: "_X" → "u_X" (ExprTk rejects '_' prefix)
    //   - Conditional:   reserved-word names registered on this evaluator
    //                    (e.g., user's `const` → `r_const`)
    std::string remapped = impl_->remap_expression(stripped);

    // Delegate to compile_preprocessed (which also caches the string)
    return compile_preprocessed(remapped);
}

int ExprTkEvaluator::compile_preprocessed(const std::string &preprocessed_expr) {
    Impl::Expression expression;
    expression.register_symbol_table(impl_->symbol_table);

    if (!impl_->parser->compile(preprocessed_expr, expression)) {
        throw std::runtime_error("ExprTk compilation failed for expression: '" + preprocessed_expr +
                                 "' — " + impl_->parser->error());
    }

    int id = static_cast<int>(impl_->expressions.size());
    impl_->expressions.push_back(std::move(expression));
    impl_->preprocessed_strings.push_back(preprocessed_expr);
    return id;
}

const std::string &ExprTkEvaluator::preprocessed_expr(int expr_id) const {
    if (expr_id < 0 || expr_id >= static_cast<int>(impl_->preprocessed_strings.size())) {
        throw std::runtime_error("Invalid expression ID: " + std::to_string(expr_id));
    }
    return impl_->preprocessed_strings[expr_id];
}

int ExprTkEvaluator::n_expressions() const { return static_cast<int>(impl_->expressions.size()); }

double ExprTkEvaluator::evaluate(int expr_id) {
    if (expr_id < 0 || expr_id >= static_cast<int>(impl_->expressions.size())) {
        throw std::runtime_error("Invalid expression ID: " + std::to_string(expr_id));
    }
    return impl_->expressions[expr_id].value();
}

void ExprTkEvaluator::set_time_ptr(double *time_addr) { impl_->set_time_ptr(time_addr); }

// ─── Efficient clone support ─────────────────────────────────────────────────

ExprTkEvaluator::ExprTkEvaluator(std::shared_ptr<void> shared_parser)
    : impl_(std::make_unique<Impl>(std::static_pointer_cast<Impl::Parser>(shared_parser))) {}

ExprTkEvaluator::ExprTkEvaluator(ExprTkEvaluator &&) noexcept = default;
ExprTkEvaluator &ExprTkEvaluator::operator=(ExprTkEvaluator &&) noexcept = default;

std::unique_ptr<ExprTkEvaluator> ExprTkEvaluator::clone_empty() const {
    // Share the parser (heavyweight, stateless between calls)
    // but create a fresh symbol table and expression list.
    return std::unique_ptr<ExprTkEvaluator>(new ExprTkEvaluator(impl_->parser));
}

// ─── Reserved names ──────────────────────────────────────────────────────────

ReservedNames reserved_names() {
    ReservedNames names;
    names.constants = {"_pi", "_e", "_kB", "_NA", "_R", "_h", "_F"};
    names.functions = {"time",  "sin",   "cos",   "tan",   "asin",  "acos", "atan",  "sinh",
                       "cosh",  "tanh",  "asinh", "acosh", "atanh", "exp",  "log",   "ln",
                       "log2",  "log10", "sqrt",  "abs",   "floor", "ceil", "round", "rint",
                       "trunc", "min",   "max",   "clamp", "avg",   "sum",  "erf",   "erfc",
                       "sign",  "sgn",   "if",    "mratio"};
    return names;
}

} // namespace bngsim
