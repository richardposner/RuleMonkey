#pragma once

#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace rulemonkey::expr {

enum class NodeType { Literal, Variable, UnaryNeg, BinaryOp, FunctionCall };

struct AstNode {
  NodeType type;
  double literal{};
  std::string name;
  char op{};
  std::vector<std::unique_ptr<AstNode>> children;

  std::unique_ptr<AstNode> clone() const;
};

using VariableMap = std::unordered_map<std::string, double>;

// Parse a BNG-style math expression into an AST.
// Throws std::runtime_error on malformed input.
std::unique_ptr<AstNode> parse(const std::string& text);

// Evaluate an AST against a variable map.
// Throws std::runtime_error if a variable is unresolved.
double evaluate(const AstNode& node, const VariableMap& vars);

// Collect all variable names referenced in the AST.
void collect_variables(const AstNode& node, std::vector<std::string>& out);

}  // namespace rulemonkey::expr
