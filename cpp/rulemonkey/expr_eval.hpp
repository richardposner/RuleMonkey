#pragma once

#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace rulemonkey::expr {

enum class NodeType { Literal, Variable, UnaryNeg, BinaryOp, FunctionCall };

// Identifies a parse-time-resolved builtin function call (math + the
// 3-arg `if` ternary).  Set by `parse()` on FunctionCall nodes whose
// name matches a known builtin.  `None` means the FunctionCall is a
// user-defined function (resolved by name lookup at evaluate time).
//
// Tagging at parse time replaces the old `is_builtin` linear scan and
// the chained `if (name == "log") ... else if (name == "ln") ...`
// dispatch in the evaluator with a single switch.
enum class BuiltinKind {
  None = 0,
  Log,
  Ln,
  Log10,
  Log2,
  Exp,
  Sqrt,
  Abs,
  Floor,
  Ceil,
  Round,
  Sin,
  Cos,
  Tan,
  Asin,
  Acos,
  Atan,
  Sinh,
  Cosh,
  Tanh,
  Min,
  Max,
  Pow,
  Atan2,
  If,
};

struct AstNode {
  NodeType type;
  double literal{};
  std::string name;
  char op{};
  std::vector<std::unique_ptr<AstNode>> children;

  // For Variable / non-builtin FunctionCall nodes: -1 until set by
  // `lower_variables()`, then the index into a flat `vector<double>`
  // representing the engine's variable layout.  When >= 0 the
  // vector-overload `evaluate(node, vars_flat)` reads `vars_flat[i]`
  // directly; `name` is kept for the legacy map evaluator and for
  // diagnostics.
  int var_index = -1;

  // For FunctionCall nodes: classified at parse time.  When != None
  // the evaluator dispatches on the enum without inspecting `name`.
  BuiltinKind builtin_kind = BuiltinKind::None;

  std::unique_ptr<AstNode> clone() const;
};

using VariableMap = std::unordered_map<std::string, double>;

// Parse a BNG-style math expression into an AST.
// Throws std::runtime_error on malformed input.
std::unique_ptr<AstNode> parse(const std::string& text);

// Evaluate an AST against a variable map.  Used by callers that don't
// have a stable variable layout (e.g., the simulator's parameter
// cascade — see resolve_cached in simulator.cpp).  Variable lookup
// goes by name through the map.
// Throws std::runtime_error if a variable is unresolved.
double evaluate(const AstNode& node, const VariableMap& vars);

// Indexed evaluator: every Variable / non-builtin FunctionCall node
// must already have `var_index >= 0` (call lower_variables() first).
// `vars_flat` must be sized to cover every var_index used by the AST.
// Used by the engine's per-event rate-law hot path; replaces a string
// hash lookup with a direct vector index.
double evaluate(const AstNode& node, const std::vector<double>& vars_flat);

// In-place AST lowering: resolve every Variable and non-builtin
// FunctionCall name to an integer index using `index_of`.  Names not
// in the map leave `var_index = -1`, and the indexed evaluator will
// throw on those.  Builtin FunctionCall nodes are left untouched (they
// dispatch on `builtin_kind`).
void lower_variables(AstNode& node, const std::unordered_map<std::string, int>& index_of);

// Collect all variable names referenced in the AST.
void collect_variables(const AstNode& node, std::vector<std::string>& out);

} // namespace rulemonkey::expr
