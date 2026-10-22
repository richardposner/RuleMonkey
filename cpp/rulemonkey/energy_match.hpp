#pragma once

// Exact free-energy changes for Arrhenius rules, counted as BioNetGen's
// network counts them (issue #109).
//
// BioNetGen weighs an energy pattern in a species by its number of matches
// divided by the pattern's automorphisms (EnergyPattern.pm, `update`), and
// gives a reaction dG = sum over patterns of G times the change in weight
// from its reactant species to its product species (`getStoich`).  For a
// connected pattern, the matches a bond changes are exactly those that put
// some pattern component on one of the bond's two components: any other
// match reads only components whose state and partner the bond leaves alone.
// So dG of one firing is the change in those anchored matches alone, counted
// on the complexes the bond joins or parts, with no species and no network.
//
// The matcher reads any graph through a small interface (a "view"), so the
// engine runs it on its molecule pool, the automorphism count runs it on a
// pattern, and the tests run it on a hand-built graph:
//
//   int  type_of(int mol) const                  molecule type index
//   int  mol_of(int comp) const                  molecule a component is on
//   int  name_of(int comp) const                 interned component name
//   int  state_of(int comp) const                internal state, -1 for none
//   int  partner_of(int comp) const              bonded component, -1 free,
//                                                kEmPhantom or kEmUnknown
//   void comps_named(int mol, int name, std::vector<int>& out) const
//                                                the molecule's components of
//                                                that name, appended to out

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace rulemonkey {

// Partner values besides a component id and -1 (free).  kEmPhantom: bound to
// something outside the graph, as a site is just before it binds a partner
// the count leaves out.  kEmUnknown: a `!?` when a pattern is the graph.
inline constexpr int kEmPhantom = -2;
inline constexpr int kEmUnknown = -3;

// A component's required state when it is `~?`: some state, any one.  A
// pattern read as a graph reports a `~?` component's state as this too.
inline constexpr int kEmAnyState = -2;

// One component of an energy pattern, resolved against its molecule type.
struct EmComponent {
  int name = -1;  // interned component name
  int state = -1; // required state index, -1 for none, kEmAnyState for `~?`
  enum class Bond { Any, Free, Bound, To };
  Bond bond = Bond::Any; // `!?` or unlisted bond, free, `!+`, explicit bond
  int to_mol = -1;       // Bond::To: the pattern component at the other end
  int to_comp = -1;
};

struct EmMolecule {
  int type = -1;
  std::vector<EmComponent> comps;
};

struct EmPattern {
  std::vector<EmMolecule> mols;
  std::string text;           // the pattern as BNGL, for messages
  double energy = 0.0;        // G, resolved
  std::string energy_expr;    // its source, re-resolved under set_param
  double automorphisms = 1.0; // matches of the pattern on itself
};

// The compiled energy patterns of a model, and the component-name table the
// engine reads molecules through.
struct ArrheniusModel {
  std::vector<EmPattern> patterns;
  std::vector<std::string> names; // interned component names
  // Each molecule type's component slots, as interned names.
  std::vector<std::vector<int>> slot_names;
  bool empty() const { return patterns.empty(); }
};

// What an Arrhenius rule that takes the exact path needs at run time (issue
// #109).  Its rate is BioNetGen's, exp(-(Ea0 + c*dG)) with c = phi (1-phi
// for the reverse half of a reversible rule) and dG counted at each firing.
// The engine samples firings against per-molecule weights that bound that
// rate from above and keeps each with probability rate / bound, so the rule
// fires at exactly the counted rate.
struct ArrheniusExact {
  bool binding = false; // makes its bond (true) or breaks it
  bool reverse_half = false;
  // The bond's two ends, as (molecule type, interned name).  For a binding
  // rule end 0 is on reactant pattern 0 and end 1 on pattern 1; for an
  // unbinding rule end 0 is on the molecule the weight is kept on.
  int type[2] = {-1, -1};
  int name[2] = {-1, -1};
  int flat[2] = {-1, -1}; // each end's flat component in the reactant pattern
  // The energy patterns the bond can change: indices into
  // ArrheniusModel::patterns.
  std::vector<int> patterns;
  // Binding only: every bond of one of those patterns that is of the rule's
  // kind, cut in two.  A match through the new bond is a match of side 0 on
  // end 0's complex times one of side 1 on end 1's (issue #109).  Where a
  // side is the cut site alone it matches exactly once, so the product is
  // the other side's matches, which that end's weight carries exactly.
  // Otherwise side 1 can match at most `bound1` ways, which is how the end-0
  // weight bounds the product.  A bond whose cut leaves the pattern
  // connected sits on a cycle, and binding two complexes cannot close one.
  struct Cut {
    int pattern = -1;
    EmPattern side0, side1;
    int root_mol0 = -1, root_comp0 = -1;   // side 0's component on end 0
    int root_mol1 = -1, root_comp1 = -1;   // side 1's component on end 1
    bool single0 = false, single1 = false; // the side is the cut site alone
    double bound1 = 0;                     // most matches side 1 can have
  };
  std::vector<Cut> cuts;
  // Binding only: the cuts with context on both sides, split by how many ways
  // side 1 matches.  This rule takes only the pairs whose end 1 matches cut
  // split[k]'s side 1 exactly split_class[k] ways; the loader makes one rule
  // per combination, so each product is exact and nothing is bounded.  A cut
  // with too many combinations to split stays bounded at bound1.
  std::vector<int> split;
  std::vector<int> split_class;
};

// ---- Views -----------------------------------------------------------------

// A view with up to two components re-partnered: a bond made or broken, or a
// site bound to a partner off the graph, without touching the graph itself.
template <class G> struct EmOverlay {
  const G& g; // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members): a view over a graph that
              // outlives it
  int a = -1, a_partner = -1;
  int b = -1, b_partner = -1;
  int type_of(int mol) const { return g.type_of(mol); }
  int mol_of(int comp) const { return g.mol_of(comp); }
  int name_of(int comp) const { return g.name_of(comp); }
  int state_of(int comp) const { return g.state_of(comp); }
  int partner_of(int comp) const {
    if (comp == a)
      return a_partner;
    if (comp == b)
      return b_partner;
    return g.partner_of(comp);
  }
  void comps_named(int mol, int name, std::vector<int>& out) const {
    g.comps_named(mol, name, out);
  }
};

// ---- Matching ----------------------------------------------------------------

namespace em_detail {

// Every match of a connected pattern that puts pattern component (pm, pc) on
// component c, each visited once.  A match maps pattern molecules to
// distinct molecules and each molecule's pattern components to distinct
// components of the same name, honouring states and bonds; it is recorded as
// the component each pattern component lands on.
template <class G> class Matcher {
public:
  Matcher(const EmPattern& p, const G& g) : p_(p), g_(g) {
    int const n = static_cast<int>(p.mols.size());
    mol_img_.assign(n, -1);
    comp_img_.resize(n);
    for (int m = 0; m < n; ++m)
      comp_img_[m].assign(p.mols[m].comps.size(), -1);
  }

  template <class F> void from(int pm, int pc, int c, F visit) {
    plan(pm);
    root_comp_ = pc;
    root_target_ = c;
    step(0, visit);
  }

  const std::vector<std::vector<int>>& images() const { return comp_img_; }
  const std::vector<int>& mol_images() const { return mol_img_; }

private:
  struct Edge {
    int parent = -1, parent_comp = -1, comp = -1;
  };

  // Molecules in the order they are matched from `root`, each after the
  // molecule whose bond reaches it.
  void plan(int root) {
    int const n = static_cast<int>(p_.mols.size());
    order_.clear();
    edge_.assign(n, Edge{});
    std::vector<char> seen(n, 0);
    order_.push_back(root);
    seen[root] = 1;
    for (std::size_t k = 0; k < order_.size(); ++k) {
      int const m = order_[k];
      const auto& comps = p_.mols[m].comps;
      for (int ci = 0; ci < static_cast<int>(comps.size()); ++ci) {
        const EmComponent& ec = comps[ci];
        if (ec.bond != EmComponent::Bond::To || seen[ec.to_mol])
          continue;
        seen[ec.to_mol] = 1;
        edge_[ec.to_mol] = Edge{m, ci, ec.to_comp};
        order_.push_back(ec.to_mol);
      }
    }
  }

  // BioNetGen's test of a pattern component on a component
  // (SpeciesGraph.pm, isomorphicToSubgraph): a required state must be there,
  // and `~?` asks only that some state is; a free site must have no bond, a
  // `!+` exactly one, which a `!+` or `!?` in a pattern read as a graph
  // counts as; `!?` asks nothing.  An explicit bond needs a real partner,
  // which the match then checks is the one the pattern names.
  bool fits(const EmComponent& ec, int c) const {
    int const state = g_.state_of(c);
    if (ec.state == kEmAnyState ? state == -1 : (ec.state >= 0 && state != ec.state))
      return false;
    int const partner = g_.partner_of(c);
    switch (ec.bond) {
    case EmComponent::Bond::Free:
      return partner == -1;
    case EmComponent::Bond::Bound:
      return partner >= 0 || partner == kEmPhantom || partner == kEmUnknown;
    case EmComponent::Bond::To:
      return partner >= 0;
    case EmComponent::Bond::Any:
      return true;
    }
    return false;
  }

  template <class F> void step(std::size_t k, F& visit) {
    if (order_.size() != p_.mols.size())
      return; // a disconnected pattern: its bonds do not reach every molecule
    if (k == order_.size()) {
      // Every pattern bond lands on a bond: the tree edges by construction,
      // the rest (a pattern with a cycle) checked here.
      for (std::size_t m = 0; m < p_.mols.size(); ++m) {
        const auto& comps = p_.mols[m].comps;
        for (std::size_t ci = 0; ci < comps.size(); ++ci)
          if (comps[ci].bond == EmComponent::Bond::To &&
              g_.partner_of(comp_img_[m][ci]) != comp_img_[comps[ci].to_mol][comps[ci].to_comp])
            return;
      }
      visit();
      return;
    }
    int const pm = order_[k];
    int mol = -1;
    int fixed_comp = -1;
    int fixed_target = -1;
    if (k == 0) {
      fixed_comp = root_comp_;
      fixed_target = root_target_;
    } else {
      const Edge& e = edge_[pm];
      int const q = g_.partner_of(comp_img_[e.parent][e.parent_comp]);
      if (q < 0)
        return;
      fixed_comp = e.comp;
      fixed_target = q;
    }
    mol = g_.mol_of(fixed_target);
    if (g_.type_of(mol) != p_.mols[pm].type ||
        g_.name_of(fixed_target) != p_.mols[pm].comps[fixed_comp].name)
      return;
    for (std::size_t j = 0; j < k; ++j)
      if (mol_img_[order_[j]] == mol)
        return;
    mol_img_[pm] = mol;
    assign(pm, 0, fixed_comp, fixed_target, k, visit);
    mol_img_[pm] = -1;
  }

  template <class F>
  void assign(int pm, int ci, int fixed_comp, int fixed_target, std::size_t k, F& visit) {
    const auto& comps = p_.mols[pm].comps;
    if (ci == static_cast<int>(comps.size())) {
      step(k + 1, visit);
      return;
    }
    const EmComponent& ec = comps[ci];
    // Candidates live in a buffer per recursion depth, read by index: deeper
    // calls may grow the outer vector and move the inner ones.
    std::size_t const depth = scratch_depth_++;
    if (depth_buf_.size() <= depth)
      depth_buf_.emplace_back();
    depth_buf_[depth].clear();
    if (ci == fixed_comp)
      depth_buf_[depth].push_back(fixed_target);
    else
      g_.comps_named(mol_img_[pm], ec.name, depth_buf_[depth]);
    for (std::size_t i = 0; i < depth_buf_[depth].size(); ++i) {
      int const c = depth_buf_[depth][i];
      bool taken = false;
      for (int j = 0; j < ci && !taken; ++j)
        taken = comp_img_[pm][j] == c;
      if (taken || !fits(ec, c))
        continue;
      comp_img_[pm][ci] = c;
      assign(pm, ci + 1, fixed_comp, fixed_target, k, visit);
      comp_img_[pm][ci] = -1;
    }
    --scratch_depth_;
  }

  // Not owned: a matcher lives for one count, inside the pattern's and the
  // graph's lifetimes.
  const EmPattern& p_; // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
  const G& g_;         // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
  std::vector<int> order_;
  std::vector<Edge> edge_;
  std::vector<int> mol_img_;
  std::vector<std::vector<int>> comp_img_;
  int root_comp_ = -1;
  int root_target_ = -1;
  std::vector<std::vector<int>> depth_buf_;
  std::size_t scratch_depth_ = 0;
};

} // namespace em_detail

// Matches of a connected pattern that put pattern component (pm, pc) on c.
template <class G> long em_count_rooted(const EmPattern& p, const G& g, int pm, int pc, int c) {
  long n = 0;
  em_detail::Matcher<G> m(p, g);
  m.from(pm, pc, c, [&n] { ++n; });
  return n;
}

// Distinct matches of a connected pattern that put some pattern component on
// one of `anchors` (component ids; negative entries are skipped).
template <class G>
long em_count_anchored(const EmPattern& p, const G& g, std::initializer_list<int> anchors) {
  std::set<std::vector<int>> seen;
  em_detail::Matcher<G> m(p, g);
  std::vector<int> key;
  for (int const a : anchors) {
    if (a < 0)
      continue;
    int const mol = g.mol_of(a);
    for (int pm = 0; pm < static_cast<int>(p.mols.size()); ++pm) {
      if (p.mols[pm].type != g.type_of(mol))
        continue;
      const auto& comps = p.mols[pm].comps;
      for (int pc = 0; pc < static_cast<int>(comps.size()); ++pc) {
        if (comps[pc].name != g.name_of(a))
          continue;
        m.from(pm, pc, a, [&] {
          key.clear();
          for (const auto& row : m.images())
            key.insert(key.end(), row.begin(), row.end());
          key.insert(key.end(), m.mol_images().begin(), m.mol_images().end());
          seen.insert(key);
        });
      }
    }
  }
  return static_cast<long>(seen.size());
}

// The free-energy change of making (make) or breaking the bond between
// components c and d: the change in BioNetGen's weights of the given
// patterns, sum of G * (matches after - matches before) / automorphisms over
// the matches anchored on c or d.
template <class G>
double em_delta_g_bond(const ArrheniusModel& am, const std::vector<int>& which, const G& g, int c,
                       int d, bool make) {
  EmOverlay<G> const after{g, c, make ? d : -1, d, make ? c : -1};
  double dg = 0.0;
  for (int const pi : which) {
    const EmPattern& p = am.patterns[pi];
    long const before = em_count_anchored(p, g, {c, d});
    long const now = em_count_anchored(p, after, {c, d});
    if (now != before)
      dg += p.energy * static_cast<double>(now - before) / p.automorphisms;
  }
  return dg;
}

// The part of a binding's free-energy change that one end's own complex
// decides: the change in the given patterns' matches anchored on component c
// when c goes from free to bound to a partner off the graph.  A match through
// the new bond itself is not among them (see ArrheniusExact::Cut).
template <class G>
double em_delta_g_end(const ArrheniusModel& am, const std::vector<int>& which, const G& g, int c) {
  EmOverlay<G> const bound{g, c, kEmPhantom};
  double dg = 0.0;
  for (int const pi : which) {
    const EmPattern& p = am.patterns[pi];
    long const before = em_count_anchored(p, g, {c});
    long const now = em_count_anchored(p, bound, {c});
    if (now != before)
      dg += p.energy * static_cast<double>(now - before) / p.automorphisms;
  }
  return dg;
}

// BioNetGen's Arrhenius rate for a firing whose reaction has free-energy
// change dg (see arrhenius_rate in energy_expand.hpp).
inline double em_arrhenius_rate(double ea0, double phi, double dg, bool reverse_half) {
  double const c = reverse_half ? 1.0 - phi : phi;
  return std::exp(-(ea0 + (c * dg)));
}

// ---- Compiling ---------------------------------------------------------------

// Whether every molecule of a pattern is reached by its bonds.
bool em_connected(const EmPattern& p);

// Whether a cut side is its cut site alone, which matches any site of its
// name exactly once.
inline bool em_single_site(const EmPattern& side) {
  return side.mols.size() == 1 && side.mols[0].comps.size() == 1 && side.mols[0].comps[0].state < 0;
}

// Matches of the pattern on itself, as BioNetGen counts its automorphisms
// (SpeciesGraph.pm: isomorphicToSubgraph of the pattern on itself).
double em_automorphisms(const EmPattern& p);

// The two sides of a pattern cut at its bond (pm, pc)-(pm2, pc2): side 0 the
// molecules reached from pm without that bond, side 1 those reached from
// pm2.  Each side's cut component keeps its state and drops its bond, since
// the bond runs to the other side.  False if the bond is on a cycle.
bool em_cut(const EmPattern& p, int pm, int pc, EmPattern& side0, int& root_mol0, int& root_comp0,
            EmPattern& side1, int& root_mol1, int& root_comp1);

// The most matches a pattern side can have with (pm, pc) on a given
// component: each molecule's image, and the component the bond reaching it
// lands on, is fixed by that bond, so only the choice among the rest of its
// same-named components multiplies.
double em_rooted_bound(const EmPattern& side, int pm, int pc, const ArrheniusModel& am);

} // namespace rulemonkey
