#include "expr_eval.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string_view>

namespace rulemonkey::expr {

// ---------------------------------------------------------------------------
// AstNode helpers
// ---------------------------------------------------------------------------

std::unique_ptr<AstNode> AstNode::clone() const {
  auto c = std::make_unique<AstNode>();
  c->type = type;
  c->literal = literal;
  c->name = name;
  c->op = op;
  c->var_index = var_index;
  c->builtin_kind = builtin_kind;
  for (auto& ch : children)
    c->children.push_back(ch->clone());
  return c;
}

// Map a function-call name to its BuiltinKind, or `None` if the name
// is not a recognized builtin.  Called once per FunctionCall node at
// parse time so the evaluator never has to compare strings.
static BuiltinKind classify_builtin(std::string_view name) {
  // 1-arg
  if (name == "log" || name == "ln")
    return BuiltinKind::Ln; // log/ln share semantics
  if (name == "log10")
    return BuiltinKind::Log10;
  if (name == "log2")
    return BuiltinKind::Log2;
  if (name == "exp")
    return BuiltinKind::Exp;
  if (name == "sqrt")
    return BuiltinKind::Sqrt;
  if (name == "abs")
    return BuiltinKind::Abs;
  if (name == "floor")
    return BuiltinKind::Floor;
  if (name == "ceil")
    return BuiltinKind::Ceil;
  if (name == "round")
    return BuiltinKind::Round;
  if (name == "sin")
    return BuiltinKind::Sin;
  if (name == "cos")
    return BuiltinKind::Cos;
  if (name == "tan")
    return BuiltinKind::Tan;
  if (name == "asin")
    return BuiltinKind::Asin;
  if (name == "acos")
    return BuiltinKind::Acos;
  if (name == "atan")
    return BuiltinKind::Atan;
  if (name == "sinh")
    return BuiltinKind::Sinh;
  if (name == "cosh")
    return BuiltinKind::Cosh;
  if (name == "tanh")
    return BuiltinKind::Tanh;
  // 2-arg
  if (name == "min")
    return BuiltinKind::Min;
  if (name == "max")
    return BuiltinKind::Max;
  if (name == "pow")
    return BuiltinKind::Pow;
  if (name == "atan2")
    return BuiltinKind::Atan2;
  // 3-arg ternary
  if (name == "if")
    return BuiltinKind::If;
  return BuiltinKind::None;
}

// ---------------------------------------------------------------------------
// Tokenizer
// ---------------------------------------------------------------------------

enum class TokType { Number, Ident, Op, LParen, RParen, Comma, End };

struct Token {
  TokType type;
  std::string text;
  double number{};
};

class Tokenizer {
public:
  explicit Tokenizer(std::string_view src) : src_(src), pos_(0) {}

  Token next() {
    skip_ws();
    if (pos_ >= src_.size())
      return {TokType::End, "", 0};

    char ch = src_[pos_];

    if (ch == '(') {
      ++pos_;
      return {TokType::LParen, "(", 0};
    }
    if (ch == ')') {
      ++pos_;
      return {TokType::RParen, ")", 0};
    }
    if (ch == ',') {
      ++pos_;
      return {TokType::Comma, ",", 0};
    }

    if (ch == '+' || ch == '-' || ch == '*' || ch == '/' || ch == '^' || ch == '<' || ch == '>' ||
        ch == '=' || ch == '!') {
      std::string op(1, ch);
      ++pos_;
      // Two-char operators: <=, >=, ==, !=, &&, ||
      if (pos_ < src_.size()) {
        char next = src_[pos_];
        if ((ch == '<' && next == '=') || (ch == '>' && next == '=') ||
            (ch == '=' && next == '=') || (ch == '!' && next == '=')) {
          op += next;
          ++pos_;
        }
      }
      return {TokType::Op, op, 0};
    }
    if (ch == '&' && pos_ + 1 < src_.size() && src_[pos_ + 1] == '&') {
      pos_ += 2;
      return {TokType::Op, "&&", 0};
    }
    if (ch == '|' && pos_ + 1 < src_.size() && src_[pos_ + 1] == '|') {
      pos_ += 2;
      return {TokType::Op, "||", 0};
    }

    if (std::isdigit(ch) || ch == '.')
      return read_number();
    if (std::isalpha(ch) || ch == '_')
      return read_ident();

    throw std::runtime_error(std::string("expr_eval: unexpected character '") + ch + "'");
  }

  Token peek() {
    auto saved = pos_;
    auto tok = next();
    pos_ = saved;
    return tok;
  }

private:
  void skip_ws() {
    while (pos_ < src_.size() && std::isspace(src_[pos_]))
      ++pos_;
  }

  Token read_number() {
    auto start = pos_;
    while (pos_ < src_.size() && (std::isdigit(src_[pos_]) || src_[pos_] == '.'))
      ++pos_;
    // Scientific notation: e/E followed by optional +/- and digits
    if (pos_ < src_.size() && (src_[pos_] == 'e' || src_[pos_] == 'E')) {
      ++pos_;
      if (pos_ < src_.size() && (src_[pos_] == '+' || src_[pos_] == '-'))
        ++pos_;
      while (pos_ < src_.size() && std::isdigit(src_[pos_]))
        ++pos_;
    }
    std::string s(src_.substr(start, pos_ - start));
    double val = std::stod(s);
    return {TokType::Number, s, val};
  }

  Token read_ident() {
    auto start = pos_;
    while (pos_ < src_.size() && (std::isalnum(src_[pos_]) || src_[pos_] == '_'))
      ++pos_;
    std::string s(src_.substr(start, pos_ - start));
    return {TokType::Ident, s, 0};
  }

  std::string_view src_;
  size_t pos_;
};

// ---------------------------------------------------------------------------
// Parser (recursive descent)
// ---------------------------------------------------------------------------

class Parser {
public:
  explicit Parser(std::string_view src) : tok_(src) { advance(); }

  std::unique_ptr<AstNode> parse_expression() {
    auto node = parse_or();
    if (cur_.type != TokType::End)
      throw std::runtime_error("expr_eval: unexpected token '" + cur_.text + "' after expression");
    return node;
  }

private:
  void advance() { cur_ = tok_.next(); }

  // or: and (|| and)*
  std::unique_ptr<AstNode> parse_or() {
    auto left = parse_and();
    while (cur_.type == TokType::Op && cur_.text == "||") {
      advance();
      auto right = parse_and();
      auto node = std::make_unique<AstNode>();
      node->type = NodeType::BinaryOp;
      node->op = '|';
      node->children.push_back(std::move(left));
      node->children.push_back(std::move(right));
      left = std::move(node);
    }
    return left;
  }

  // and: comparison (&& comparison)*
  std::unique_ptr<AstNode> parse_and() {
    auto left = parse_comparison();
    while (cur_.type == TokType::Op && cur_.text == "&&") {
      advance();
      auto right = parse_comparison();
      auto node = std::make_unique<AstNode>();
      node->type = NodeType::BinaryOp;
      node->op = '&';
      node->children.push_back(std::move(left));
      node->children.push_back(std::move(right));
      left = std::move(node);
    }
    return left;
  }

  // comparison: add_sub ((<|>|<=|>=|==|!=) add_sub)?
  std::unique_ptr<AstNode> parse_comparison() {
    auto left = parse_add_sub();
    if (cur_.type == TokType::Op && (cur_.text == "<" || cur_.text == ">" || cur_.text == "<=" ||
                                     cur_.text == ">=" || cur_.text == "==" || cur_.text == "!=")) {
      char op;
      if (cur_.text == "<")
        op = '<';
      else if (cur_.text == ">")
        op = '>';
      else if (cur_.text == "<=")
        op = 'l'; // l for less-equal
      else if (cur_.text == ">=")
        op = 'g'; // g for greater-equal
      else if (cur_.text == "==")
        op = '=';
      else
        op = '!'; // !=
      advance();
      auto right = parse_add_sub();
      auto node = std::make_unique<AstNode>();
      node->type = NodeType::BinaryOp;
      node->op = op;
      node->children.push_back(std::move(left));
      node->children.push_back(std::move(right));
      return node;
    }
    return left;
  }

  // add_sub: mul_div ((+|-) mul_div)*
  std::unique_ptr<AstNode> parse_add_sub() {
    auto left = parse_mul_div();
    while (cur_.type == TokType::Op && (cur_.text == "+" || cur_.text == "-")) {
      char op = cur_.text[0];
      advance();
      auto right = parse_mul_div();
      auto node = std::make_unique<AstNode>();
      node->type = NodeType::BinaryOp;
      node->op = op;
      node->children.push_back(std::move(left));
      node->children.push_back(std::move(right));
      left = std::move(node);
    }
    return left;
  }

  // mul_div: power ((*|/) power)*
  std::unique_ptr<AstNode> parse_mul_div() {
    auto left = parse_power();
    while (cur_.type == TokType::Op && (cur_.text == "*" || cur_.text == "/")) {
      char op = cur_.text[0];
      advance();
      auto right = parse_power();
      auto node = std::make_unique<AstNode>();
      node->type = NodeType::BinaryOp;
      node->op = op;
      node->children.push_back(std::move(left));
      node->children.push_back(std::move(right));
      left = std::move(node);
    }
    return left;
  }

  // power: unary (^ power)?   (right-associative)
  std::unique_ptr<AstNode> parse_power() {
    auto left = parse_unary();
    if (cur_.type == TokType::Op && cur_.text == "^") {
      advance();
      auto right = parse_power(); // right-recursive for right-assoc
      auto node = std::make_unique<AstNode>();
      node->type = NodeType::BinaryOp;
      node->op = '^';
      node->children.push_back(std::move(left));
      node->children.push_back(std::move(right));
      return node;
    }
    return left;
  }

  // unary: (-|+)? power-operand
  //
  // The `-` / `+` operand is a full power-expression (parse_power) — not a
  // bare primary — so unary binds LOOSER than `^`.  That makes `-2^2` parse
  // as `-(2^2) = -4`, matching Python / BNG2 / standard math.  Reading the
  // primary directly here would produce `(-2)^2 = 4`, a silent wrong answer
  // for any rate law or TFUN expression with a negative base.
  //
  // The fallthrough (no leading sign) still calls parse_primary directly,
  // which is what parse_power's LHS expects — wiring it through parse_power
  // here would infinite-recurse since parse_power's LHS is parse_unary.
  std::unique_ptr<AstNode> parse_unary() {
    if (cur_.type == TokType::Op && cur_.text == "-") {
      advance();
      auto child = parse_power();
      auto node = std::make_unique<AstNode>();
      node->type = NodeType::UnaryNeg;
      node->children.push_back(std::move(child));
      return node;
    }
    if (cur_.type == TokType::Op && cur_.text == "+") {
      advance();
      return parse_power();
    }
    return parse_primary();
  }

  // primary: number | ident | ident(args) | (expr)
  std::unique_ptr<AstNode> parse_primary() {
    if (cur_.type == TokType::Number) {
      auto node = std::make_unique<AstNode>();
      node->type = NodeType::Literal;
      node->literal = cur_.number;
      advance();
      return node;
    }

    if (cur_.type == TokType::Ident) {
      std::string name = cur_.text;
      advance();

      // Check for function call: ident(...)
      if (cur_.type == TokType::LParen) {
        advance(); // consume '('
        auto node = std::make_unique<AstNode>();
        node->type = NodeType::FunctionCall;
        node->name = name;
        node->builtin_kind = classify_builtin(name);

        if (cur_.type != TokType::RParen) {
          node->children.push_back(parse_or());
          while (cur_.type == TokType::Comma) {
            advance();
            node->children.push_back(parse_or());
          }
        }
        if (cur_.type != TokType::RParen)
          throw std::runtime_error("expr_eval: expected ')' after function arguments");
        advance();
        return node;
      }

      // Named constants
      if (name == "pi" || name == "PI") {
        auto node = std::make_unique<AstNode>();
        node->type = NodeType::Literal;
        node->literal = M_PI;
        return node;
      }
      if (name == "e" || name == "E") {
        auto node = std::make_unique<AstNode>();
        node->type = NodeType::Literal;
        node->literal = M_E;
        return node;
      }

      // Variable reference
      auto node = std::make_unique<AstNode>();
      node->type = NodeType::Variable;
      node->name = name;
      return node;
    }

    if (cur_.type == TokType::LParen) {
      advance();
      auto node = parse_or();
      if (cur_.type != TokType::RParen)
        throw std::runtime_error("expr_eval: expected ')'");
      advance();
      return node;
    }

    throw std::runtime_error("expr_eval: unexpected token '" + cur_.text + "'");
  }

  Tokenizer tok_;
  Token cur_;
};

// ---------------------------------------------------------------------------
// Public parse
// ---------------------------------------------------------------------------

std::unique_ptr<AstNode> parse(const std::string& text) {
  if (text.empty())
    throw std::runtime_error("expr_eval: empty expression");
  Parser p(text);
  return p.parse_expression();
}

// ---------------------------------------------------------------------------
// Evaluator
// ---------------------------------------------------------------------------

// Builtin dispatch keyed on `BuiltinKind` — set once at parse time by
// `classify_builtin`.  Returns std::nullopt when the kind has no
// signature for the given arity (which is a user error: `if(x, y)`
// was almost certainly meant as the 3-arg ternary, not as `if`
// shadowed by a global function).  Result may itself be NaN — that's
// the well-defined IEEE-754 result for out-of-domain input
// (pow(-2, 0.5), acos(2), log(-1), …) and is propagated as-is.
//
// Templated on `EvalChild` so both the map evaluator (legacy, used by
// the simulator parameter cascade) and the indexed evaluator (engine
// per-event hot path) share the same dispatch chain.
template <class EvalChild>
static std::optional<double> eval_builtin(BuiltinKind kind,
                                          const std::vector<std::unique_ptr<AstNode>>& args,
                                          EvalChild eval_child) {
  auto n = args.size();

  // 3-arg ternary: lazy-evaluate — the branch not taken must not run
  // (e.g., `if(x>0, 1/x, 0)` must not divide by zero when x==0).
  if (kind == BuiltinKind::If && n == 3) {
    double cond = eval_child(*args[0]);
    return (cond != 0.0) ? eval_child(*args[1]) : eval_child(*args[2]);
  }

  if (n == 1) {
    double a = eval_child(*args[0]);
    switch (kind) {
    case BuiltinKind::Ln:
      return std::log(a);
    case BuiltinKind::Log10:
      return std::log10(a);
    case BuiltinKind::Log2:
      return std::log2(a);
    case BuiltinKind::Exp:
      return std::exp(a);
    case BuiltinKind::Sqrt:
      return std::sqrt(a);
    case BuiltinKind::Abs:
      return std::fabs(a);
    case BuiltinKind::Floor:
      return std::floor(a);
    case BuiltinKind::Ceil:
      return std::ceil(a);
    case BuiltinKind::Round:
      return std::round(a);
    case BuiltinKind::Sin:
      return std::sin(a);
    case BuiltinKind::Cos:
      return std::cos(a);
    case BuiltinKind::Tan:
      return std::tan(a);
    case BuiltinKind::Asin:
      return std::asin(a);
    case BuiltinKind::Acos:
      return std::acos(a);
    case BuiltinKind::Atan:
      return std::atan(a);
    case BuiltinKind::Sinh:
      return std::sinh(a);
    case BuiltinKind::Cosh:
      return std::cosh(a);
    case BuiltinKind::Tanh:
      return std::tanh(a);
    default:
      break;
    }
  }

  if (n == 2) {
    double a = eval_child(*args[0]);
    double b = eval_child(*args[1]);
    switch (kind) {
    case BuiltinKind::Min:
      return std::fmin(a, b);
    case BuiltinKind::Max:
      return std::fmax(a, b);
    case BuiltinKind::Pow:
      return std::pow(a, b);
    case BuiltinKind::Atan2:
      return std::atan2(a, b);
    default:
      break;
    }
  }

  return std::nullopt;
}

// Templated evaluator core.  `lookup_var(node)` is the only piece
// that differs between the map and indexed overloads — for Variable
// nodes it returns the variable's current value (or throws on
// unresolved); for non-builtin FunctionCall nodes we treat `f()` as
// the value of the variable `f` (this is how BNG global functions
// are referenced — the caller arranges the var lookup to return the
// function's evaluated value).
template <class LookupFn> static double evaluate_impl(const AstNode& node, LookupFn lookup_var) {
  switch (node.type) {
  case NodeType::Literal:
    return node.literal;

  case NodeType::Variable:
    return lookup_var(node);

  case NodeType::UnaryNeg:
    return -evaluate_impl(*node.children[0], lookup_var);

  case NodeType::BinaryOp: {
    // Short-circuit `&&` and `||` BEFORE evaluating the RHS —
    // `if(x>0 && 1/x>1, ...)` must not divide by zero when x==0,
    // and a NaN injected from the unreachable branch silently
    // poisons rate-law evaluation downstream.  Comparison ops
    // produce 0.0/1.0; we coerce the RHS the same way so the
    // result is always 0.0/1.0 regardless of which side decided.
    if (node.op == '&') {
      double l = evaluate_impl(*node.children[0], lookup_var);
      if (l == 0.0)
        return 0.0;
      double r = evaluate_impl(*node.children[1], lookup_var);
      return (r != 0.0) ? 1.0 : 0.0;
    }
    if (node.op == '|') {
      double l = evaluate_impl(*node.children[0], lookup_var);
      if (l != 0.0)
        return 1.0;
      double r = evaluate_impl(*node.children[1], lookup_var);
      return (r != 0.0) ? 1.0 : 0.0;
    }
    double l = evaluate_impl(*node.children[0], lookup_var);
    double r = evaluate_impl(*node.children[1], lookup_var);
    switch (node.op) {
    case '+':
      return l + r;
    case '-':
      return l - r;
    case '*':
      return l * r;
    case '/':
      return l / r;
    case '^':
      return std::pow(l, r);
    case '<':
      return (l < r) ? 1.0 : 0.0;
    case '>':
      return (l > r) ? 1.0 : 0.0;
    case 'l':
      return (l <= r) ? 1.0 : 0.0; // <=
    case 'g':
      return (l >= r) ? 1.0 : 0.0; // >=
    // == and != on doubles use bit-exact comparison, matching NFsim's
    // muParser-backed FuncFactory (which exposes the muParser default
    // operator==).  Two consequences callers must understand:
    //   (a) `if(x == 1.0, …)` after FP arithmetic that "should" land on
    //       1.0 (e.g. `x = 0.1 + 0.2`) will compare false because the
    //       result is 0.30000000000000004, not 0.3.  Authors must avoid
    //       == on derived values; use tolerance windows in BNGL itself
    //       (`if(x > 0.99 && x < 1.01, …)`) where exact equality cannot
    //       be guaranteed.
    //   (b) Two identical literal expressions evaluated via the same
    //       AST will compare equal because they produce bit-identical
    //       double results — the parse-once / lower-once pipeline does
    //       not introduce path-dependent rounding.
    // We deliberately do NOT switch to `std::abs(l-r) <= eps * std::max(...)`
    // — that is a behavior change versus NFsim and would silently shift
    // results for any model that depends on integer counts compared
    // exactly (a common pattern in hand-authored BNGL).
    case '=':
      return (l == r) ? 1.0 : 0.0; // ==
    case '!':
      return (l != r) ? 1.0 : 0.0; // !=
    default:
      throw std::runtime_error(std::string("expr_eval: unknown binary op '") + node.op + "'");
    }
  }

  case NodeType::FunctionCall: {
    if (node.builtin_kind != BuiltinKind::None) {
      auto result =
          eval_builtin(node.builtin_kind, node.children, [&lookup_var](const AstNode& c) -> double {
            return evaluate_impl(c, lookup_var);
          });
      if (result)
        return *result;
      throw std::runtime_error("expr_eval: builtin '" + node.name + "' called with " +
                               std::to_string(node.children.size()) +
                               " argument(s); no matching signature");
    }
    // Non-builtin: BNG global function reference.
    return lookup_var(node);
  }
  }
  throw std::runtime_error("expr_eval: unknown node type");
}

double evaluate(const AstNode& node, const VariableMap& vars) {
  return evaluate_impl(node, [&vars](const AstNode& n) -> double {
    auto it = vars.find(n.name);
    if (it == vars.end()) {
      if (n.type == NodeType::FunctionCall) {
        throw std::runtime_error("expr_eval: unknown function '" + n.name + "' with " +
                                 std::to_string(n.children.size()) + " arguments");
      }
      throw std::runtime_error("expr_eval: unresolved variable '" + n.name + "'");
    }
    return it->second;
  });
}

double evaluate(const AstNode& node, const std::vector<double>& vars_flat) {
  return evaluate_impl(node, [&vars_flat](const AstNode& n) -> double {
    int idx = n.var_index;
    if (idx < 0 || static_cast<size_t>(idx) >= vars_flat.size()) {
      if (n.type == NodeType::FunctionCall) {
        throw std::runtime_error("expr_eval: unknown function '" + n.name + "' with " +
                                 std::to_string(n.children.size()) + " arguments");
      }
      throw std::runtime_error("expr_eval: unresolved variable '" + n.name + "'");
    }
    return vars_flat[idx];
  });
}

void lower_variables(AstNode& node, const std::unordered_map<std::string, int>& index_of) {
  if (node.type == NodeType::Variable ||
      (node.type == NodeType::FunctionCall && node.builtin_kind == BuiltinKind::None)) {
    auto it = index_of.find(node.name);
    node.var_index = (it != index_of.end()) ? it->second : -1;
  }
  for (auto& ch : node.children)
    lower_variables(*ch, index_of);
}

// ---------------------------------------------------------------------------
// Variable collection
// ---------------------------------------------------------------------------

void collect_variables(const AstNode& node, std::vector<std::string>& out) {
  if (node.type == NodeType::Variable) {
    if (std::find(out.begin(), out.end(), node.name) == out.end())
      out.push_back(node.name);
    return;
  }
  if (node.type == NodeType::FunctionCall && node.builtin_kind == BuiltinKind::None) {
    // Non-builtin function call — the function name is a dependency
    if (std::find(out.begin(), out.end(), node.name) == out.end())
      out.push_back(node.name);
  }
  for (auto& ch : node.children)
    collect_variables(*ch, out);
}

} // namespace rulemonkey::expr
