#include "expr_eval.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
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
  for (auto& ch : children) c->children.push_back(ch->clone());
  return c;
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
    if (pos_ >= src_.size()) return {TokType::End, "", 0};

    char ch = src_[pos_];

    if (ch == '(') { ++pos_; return {TokType::LParen, "(", 0}; }
    if (ch == ')') { ++pos_; return {TokType::RParen, ")", 0}; }
    if (ch == ',') { ++pos_; return {TokType::Comma, ",", 0}; }

    if (ch == '+' || ch == '-' || ch == '*' || ch == '/' || ch == '^' ||
        ch == '<' || ch == '>' || ch == '=' || ch == '!') {
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

    if (std::isdigit(ch) || ch == '.') return read_number();
    if (std::isalpha(ch) || ch == '_') return read_ident();

    throw std::runtime_error(std::string("expr_eval: unexpected character '") +
                             ch + "'");
  }

  Token peek() {
    auto saved = pos_;
    auto tok = next();
    pos_ = saved;
    return tok;
  }

 private:
  void skip_ws() {
    while (pos_ < src_.size() && std::isspace(src_[pos_])) ++pos_;
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
      while (pos_ < src_.size() && std::isdigit(src_[pos_])) ++pos_;
    }
    std::string s(src_.substr(start, pos_ - start));
    double val = std::stod(s);
    return {TokType::Number, s, val};
  }

  Token read_ident() {
    auto start = pos_;
    while (pos_ < src_.size() &&
           (std::isalnum(src_[pos_]) || src_[pos_] == '_'))
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
      throw std::runtime_error("expr_eval: unexpected token '" + cur_.text +
                               "' after expression");
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
    if (cur_.type == TokType::Op &&
        (cur_.text == "<" || cur_.text == ">" || cur_.text == "<=" ||
         cur_.text == ">=" || cur_.text == "==" || cur_.text == "!=")) {
      char op;
      if (cur_.text == "<") op = '<';
      else if (cur_.text == ">") op = '>';
      else if (cur_.text == "<=") op = 'l';  // l for less-equal
      else if (cur_.text == ">=") op = 'g';  // g for greater-equal
      else if (cur_.text == "==") op = '=';
      else op = '!';  // !=
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
    while (cur_.type == TokType::Op &&
           (cur_.text == "+" || cur_.text == "-")) {
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
    while (cur_.type == TokType::Op &&
           (cur_.text == "*" || cur_.text == "/")) {
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
      auto right = parse_power();  // right-recursive for right-assoc
      auto node = std::make_unique<AstNode>();
      node->type = NodeType::BinaryOp;
      node->op = '^';
      node->children.push_back(std::move(left));
      node->children.push_back(std::move(right));
      return node;
    }
    return left;
  }

  // unary: (-|+)? primary
  std::unique_ptr<AstNode> parse_unary() {
    if (cur_.type == TokType::Op && cur_.text == "-") {
      advance();
      auto child = parse_primary();
      auto node = std::make_unique<AstNode>();
      node->type = NodeType::UnaryNeg;
      node->children.push_back(std::move(child));
      return node;
    }
    if (cur_.type == TokType::Op && cur_.text == "+") {
      advance();
      return parse_primary();
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
        advance();  // consume '('
        auto node = std::make_unique<AstNode>();
        node->type = NodeType::FunctionCall;
        node->name = name;

        if (cur_.type != TokType::RParen) {
          node->children.push_back(parse_or());
          while (cur_.type == TokType::Comma) {
            advance();
            node->children.push_back(parse_or());
          }
        }
        if (cur_.type != TokType::RParen)
          throw std::runtime_error(
              "expr_eval: expected ')' after function arguments");
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

static double eval_builtin(const std::string& name,
                           const std::vector<std::unique_ptr<AstNode>>& args,
                           const VariableMap& vars) {
  auto n = args.size();

  // 1-arg functions
  if (n == 1) {
    double a = evaluate(*args[0], vars);
    if (name == "log" || name == "ln") return std::log(a);
    if (name == "log10") return std::log10(a);
    if (name == "log2") return std::log2(a);
    if (name == "exp") return std::exp(a);
    if (name == "sqrt") return std::sqrt(a);
    if (name == "abs") return std::fabs(a);
    if (name == "floor") return std::floor(a);
    if (name == "ceil") return std::ceil(a);
    if (name == "round") return std::round(a);
    if (name == "sin") return std::sin(a);
    if (name == "cos") return std::cos(a);
    if (name == "tan") return std::tan(a);
    if (name == "asin") return std::asin(a);
    if (name == "acos") return std::acos(a);
    if (name == "atan") return std::atan(a);
    if (name == "sinh") return std::sinh(a);
    if (name == "cosh") return std::cosh(a);
    if (name == "tanh") return std::tanh(a);
  }

  // 2-arg functions
  if (n == 2) {
    double a = evaluate(*args[0], vars);
    double b = evaluate(*args[1], vars);
    if (name == "min") return std::fmin(a, b);
    if (name == "max") return std::fmax(a, b);
    if (name == "pow") return std::pow(a, b);
    if (name == "atan2") return std::atan2(a, b);
  }

  // if(cond, then, else) — lazy evaluation
  if (name == "if" && n == 3) {
    double cond = evaluate(*args[0], vars);
    return (cond != 0.0) ? evaluate(*args[1], vars)
                         : evaluate(*args[2], vars);
  }

  return std::numeric_limits<double>::quiet_NaN();  // sentinel: not a builtin
}

static bool is_builtin(const std::string& name) {
  static const std::vector<std::string> names = {
      "log",  "ln",   "log10", "log2", "exp",  "sqrt", "abs",
      "floor", "ceil", "round", "sin",  "cos",  "tan",  "asin",
      "acos", "atan", "sinh",  "cosh", "tanh", "min",  "max",
      "pow",  "atan2", "if"};
  return std::find(names.begin(), names.end(), name) != names.end();
}

double evaluate(const AstNode& node, const VariableMap& vars) {
  switch (node.type) {
    case NodeType::Literal:
      return node.literal;

    case NodeType::Variable: {
      auto it = vars.find(node.name);
      if (it == vars.end())
        throw std::runtime_error("expr_eval: unresolved variable '" +
                                 node.name + "'");
      return it->second;
    }

    case NodeType::UnaryNeg:
      return -evaluate(*node.children[0], vars);

    case NodeType::BinaryOp: {
      double l = evaluate(*node.children[0], vars);
      double r = evaluate(*node.children[1], vars);
      switch (node.op) {
        case '+': return l + r;
        case '-': return l - r;
        case '*': return l * r;
        case '/': return l / r;
        case '^': return std::pow(l, r);
        case '<': return (l < r) ? 1.0 : 0.0;
        case '>': return (l > r) ? 1.0 : 0.0;
        case 'l': return (l <= r) ? 1.0 : 0.0;  // <=
        case 'g': return (l >= r) ? 1.0 : 0.0;  // >=
        case '=': return (l == r) ? 1.0 : 0.0;  // ==
        case '!': return (l != r) ? 1.0 : 0.0;  // !=
        case '&': return (l != 0.0 && r != 0.0) ? 1.0 : 0.0;  // &&
        case '|': return (l != 0.0 || r != 0.0) ? 1.0 : 0.0;  // ||
        default:
          throw std::runtime_error(
              std::string("expr_eval: unknown binary op '") + node.op + "'");
      }
    }

    case NodeType::FunctionCall: {
      // Try builtin first
      if (is_builtin(node.name)) {
        double result = eval_builtin(node.name, node.children, vars);
        if (!std::isnan(result) || node.name == "log" || node.name == "ln" ||
            node.name == "sqrt") {
          // NaN is a valid result from log(-1) etc.
          return result;
        }
      }

      // Fall through to variable-as-function: treat f() as variable "f"
      // This supports BNG global functions referenced in expressions
      auto it = vars.find(node.name);
      if (it != vars.end()) return it->second;

      throw std::runtime_error("expr_eval: unknown function '" + node.name +
                               "' with " + std::to_string(node.children.size()) +
                               " arguments");
    }
  }
  throw std::runtime_error("expr_eval: unknown node type");
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
  if (node.type == NodeType::FunctionCall && !is_builtin(node.name)) {
    // Non-builtin function call — the function name is a dependency
    if (std::find(out.begin(), out.end(), node.name) == out.end())
      out.push_back(node.name);
  }
  for (auto& ch : node.children) collect_variables(*ch, out);
}

}  // namespace rulemonkey::expr
