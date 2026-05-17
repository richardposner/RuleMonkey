// Canonical complex labeling — implementation of the pure core.
//
// See canonical.hpp and dev/canonical_labeling_plan_2026_05_16.md.
//
// Algorithm (plan §3):
//
//   Encoding (§3.1) — a vertex-colored port graph:
//     * one vertex per molecule; its color folds in the molecule type
//       name and the sorted multiset of UNBONDED (name, state) pairs;
//     * one vertex per BONDED component; its color is (name, state);
//     * edges: molecule <-> each of its bonded-component vertices, and
//       one edge per bond between the two bonded-component vertices.
//     Folding unbonded components into the molecule color (rather than
//     giving every component a vertex, as NFsim does) keeps the graph
//     small: unbonded components carry no edges, so they cannot affect
//     graph structure.
//
//   Canonicalization (§3.2):
//     1. 1-WL color refinement — recolor each vertex by (own color,
//        sorted neighbor colors) until the partition stabilizes.
//     2. Fast path — when refinement gives every molecule vertex its
//        own color AND no molecule has two interchangeable bonded
//        components left tied, the refined colors fix a unique,
//        isomorphism-invariant ordering; render directly, no search.
//        This is the overwhelmingly common case (asymmetric complexes).
//     3. Individualization-refinement (plan §3.2 step 3) for residual
//        symmetric classes — rings, homo-oligomers.  Pick a target
//        cell, individualize each vertex in it in turn, re-refine,
//        recurse to a leaf (a coloring the fast-path test accepts);
//        the canonical label is the lexicographically minimal rendered
//        BNGL string over all leaves.  This makes the label a *true*
//        canonical form on every input — 1-WL alone is incomplete on
//        regular/symmetric graphs.
//
//   Rendering (§3.3): walk molecules in canonical order; within each
//   molecule keep the molecule type's component-name layout but order
//   interchangeable same-name components by their refined color; emit
//   BNGL `Type(comp~state!bond,...)` with bond labels assigned in
//   first-encounter order along the walk.  The string is both the dedup
//   key and the `.species` line.
//
//   Component-order contract: the renderer preserves the *name layout*
//   of the input ComplexGraph (slot i keeps its declared name) and only
//   reorders components that share a name.  `extract_complex` must
//   therefore emit each molecule's components in molecule-type
//   declaration order — which it does, since the engine stores
//   MoleculeInstance::comp_ids that way.

#include "canonical.hpp"

#include <algorithm>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace rulemonkey::canonical {

// ===========================================================================
// ComplexGraph builder
// ===========================================================================

int ComplexGraph::add_molecule(const std::string& type_name,
                               const std::vector<std::pair<std::string, std::string>>& comps) {
  Molecule m;
  m.type_name = type_name;
  m.first_comp = static_cast<int>(components.size());
  m.n_comp = static_cast<int>(comps.size());
  for (const auto& [name, state] : comps)
    components.push_back(Component{name, state, -1});
  molecules.push_back(std::move(m));
  return static_cast<int>(molecules.size()) - 1;
}

void ComplexGraph::add_bond(int mol_a, int local_a, int mol_b, int local_b) {
  int const ga = global_comp(mol_a, local_a);
  int const gb = global_comp(mol_b, local_b);
  components[ga].partner = gb;
  components[gb].partner = ga;
}

// ===========================================================================
// Canonicalization
// ===========================================================================

namespace {

// Leaf-count safety valve for individualization-refinement (plan §3.2
// step 3, and the step-3 brief's "recursion/leaf guard").  BNGL
// complexes are small and richly labeled, so the search tree is tiny
// and shallow in practice; this cap exists only so a pathological
// hand-built input cannot blow up.  It is NOT a tuning knob — plan §5
// defers all profiling/tuning to step 7.  If the cap is ever hit the
// best leaf found so far is returned (still deterministic per input).
constexpr long kLeafBudget = 200000;

// A bonded component contributes a vertex whose initial color is this
// string.  Unbonded components instead fold into their molecule's color
// (see molecule_color).
std::string component_color(const ComplexGraph::Component& c) {
  return "C:" + c.name + "~" + c.state;
}

// A molecule's initial color: the type name plus the sorted multiset of
// its UNBONDED components' `(name, state)`.  Bonded components are
// excluded — they become their own vertices and carry the edges.  The
// `M:` / `C:` prefixes keep molecule and component colors in disjoint
// ranges of the interning table.
std::string molecule_color(const ComplexGraph& g, const ComplexGraph::Molecule& m) {
  std::vector<std::string> free_comps;
  for (int i = 0; i < m.n_comp; ++i) {
    const auto& c = g.components[m.first_comp + i];
    if (c.partner < 0)
      free_comps.push_back(c.name + "~" + c.state);
  }
  std::sort(free_comps.begin(), free_comps.end());
  std::string s = "M:" + m.type_name + "(";
  for (size_t i = 0; i < free_comps.size(); ++i) {
    if (i)
      s += ',';
    s += free_comps[i];
  }
  s += ')';
  return s;
}

// Within-molecule canonical ordering key for one component.  Components
// that share a name are interchangeable BNGL sites; the renderer orders
// each such group by this key.  Free components (which have no graph
// vertex) sort ahead of bonded ones and break ties by state; bonded
// components sort by their refined WL color.  `bonded_color` is the
// refined color of the component's vertex, or -1 for a free component.
using CompKey = std::tuple<int, int, std::string>; // (bonded?, color, state)
CompKey comp_key(const ComplexGraph::Component& c, int bonded_color) {
  if (bonded_color < 0)
    return {0, 0, c.state};
  return {1, bonded_color, std::string{}};
}

// 1-WL color refinement, in place.  Each round recolors a vertex by
// (own color, sorted neighbor colors); signatures are sorted and
// ranked, so the output colors are a canonical 0..k-1 ranking derived
// purely from structure.  The own-color term means refinement only
// ever splits classes — it has converged once the class count stops
// rising.  Accepts any integer colors on input (individualization
// passes a vertex an out-of-range marker color; the first round
// re-ranks everything), so it doubles as the post-individualization
// re-refine step.
void refine(const std::vector<std::vector<int>>& adj, std::vector<int>& color) {
  int const n_vert = static_cast<int>(color.size());
  if (n_vert == 0)
    return;
  std::vector<int> next_color(n_vert);
  int n_classes = 0;
  bool first = true;
  while (true) {
    std::vector<std::pair<std::vector<int>, int>> sigs; // (signature, vertex)
    sigs.reserve(n_vert);
    for (int v = 0; v < n_vert; ++v) {
      std::vector<int> sig;
      sig.reserve(adj[v].size() + 1);
      sig.push_back(color[v]);
      std::vector<int> nbr;
      nbr.reserve(adj[v].size());
      for (int const u : adj[v])
        nbr.push_back(color[u]);
      std::sort(nbr.begin(), nbr.end());
      sig.insert(sig.end(), nbr.begin(), nbr.end());
      sigs.emplace_back(std::move(sig), v);
    }
    std::sort(sigs.begin(), sigs.end());

    int rank = -1;
    const std::vector<int>* prev = nullptr;
    for (const auto& [sig, v] : sigs) {
      if (prev == nullptr || *prev != sig)
        ++rank;
      next_color[v] = rank;
      prev = &sig;
    }
    int const new_classes = rank + 1;
    color.swap(next_color);
    if (!first && new_classes == n_classes)
      break; // partition stable
    first = false;
    n_classes = new_classes;
  }
}

// Refined color of each bonded component's vertex (-1 if free) — the
// renderer's within-molecule ordering key.
std::vector<int> bonded_colors(const ComplexGraph& g, const std::vector<int>& color,
                               const std::vector<int>& comp_vertex) {
  std::vector<int> bc(g.components.size(), -1);
  for (int gc = 0; gc < static_cast<int>(g.components.size()); ++gc) {
    if (comp_vertex[gc] >= 0)
      bc[gc] = color[comp_vertex[gc]];
  }
  return bc;
}

// Leaf test: is `color` discriminating enough that `render` produces a
// unique, isomorphism-invariant string?  True iff (a) every molecule
// vertex landed in its own color class and (b) no molecule has two
// interchangeable bonded components — same name, same refined color —
// left tied.  Either residual symmetry is what individualization
// resolves; this same test is the fast-path gate and the search's leaf
// condition.  (Component vertices on *different* molecules, or of
// different names, may still share a color — that never affects the
// render, so it does not block a leaf.)
bool is_leaf(const ComplexGraph& g, const std::vector<int>& color,
             const std::vector<int>& comp_vertex, int n_mol) {
  std::vector<int> mc(color.begin(), color.begin() + n_mol);
  std::sort(mc.begin(), mc.end());
  for (int i = 1; i < n_mol; ++i)
    if (mc[i] == mc[i - 1])
      return false;
  for (int m = 0; m < n_mol; ++m) {
    const auto& mol = g.molecules[m];
    std::map<std::string, std::vector<int>> seen; // name -> bonded-vertex colors
    for (int i = 0; i < mol.n_comp; ++i) {
      int const gc = mol.first_comp + i;
      if (comp_vertex[gc] < 0)
        continue; // free component — no vertex, physically interchangeable
      int const c = color[comp_vertex[gc]];
      auto& bucket = seen[g.components[gc].name];
      for (int const pc : bucket)
        if (pc == c)
          return false;
      bucket.push_back(c);
    }
  }
  return true;
}

// Render the canonical BNGL string.  `order` is the canonical molecule
// ordering; `bonded_color[gc]` is the refined color of bonded component
// `gc`'s vertex (-1 if `gc` is free).  Within each molecule the
// type's name layout is preserved and same-name components are placed
// in `comp_key` order; bond labels are assigned 1,2,3,... in
// first-encounter order along the resulting walk.
std::string render(const ComplexGraph& g, const std::vector<int>& order,
                   const std::vector<int>& bonded_color) {
  std::vector<int> comp_label(g.components.size(), 0); // 0 = unassigned
  int next_label = 1;
  std::string out;

  for (size_t oi = 0; oi < order.size(); ++oi) {
    if (oi)
      out += '.';
    const auto& m = g.molecules[order[oi]];

    // Group this molecule's component slots by name, then order each
    // group canonically.  `groups[name]` is the global component
    // indices of that name, sorted by comp_key.
    std::map<std::string, std::vector<int>> groups;
    for (int i = 0; i < m.n_comp; ++i)
      groups[g.components[m.first_comp + i].name].push_back(m.first_comp + i);
    for (auto& [name, gcs] : groups) {
      std::sort(gcs.begin(), gcs.end(), [&](int a, int b) {
        return comp_key(g.components[a], bonded_color[a]) <
               comp_key(g.components[b], bonded_color[b]);
      });
    }
    std::map<std::string, int> cursor;

    out += m.type_name;
    out += '(';
    for (int i = 0; i < m.n_comp; ++i) {
      if (i)
        out += ',';
      // Slot i keeps its declared name; the canonical component placed
      // there is the next one from that name's ordered group.
      const std::string& name = g.components[m.first_comp + i].name;
      int const gc = groups[name][cursor[name]++];
      const auto& c = g.components[gc];
      out += c.name;
      if (!c.state.empty()) {
        out += '~';
        out += c.state;
      }
      if (c.partner >= 0) {
        int& lbl = comp_label[gc];
        if (lbl == 0) {
          lbl = next_label++;
          comp_label[c.partner] = lbl;
        }
        out += '!';
        out += std::to_string(lbl);
      }
    }
    out += ')';
  }
  return out;
}

// Render a leaf coloring: derive the canonical molecule order (sort by
// refined color — at a leaf every molecule color is distinct, so this
// is a total order) and the per-component bonded colors, then render.
std::string render_leaf(const ComplexGraph& g, const std::vector<int>& color,
                        const std::vector<int>& comp_vertex, int n_mol) {
  std::vector<int> order(n_mol);
  for (int m = 0; m < n_mol; ++m)
    order[m] = m;
  std::sort(order.begin(), order.end(), [&](int a, int b) { return color[a] < color[b]; });
  return render(g, order, bonded_colors(g, color, comp_vertex));
}

// Pick the target cell for individualization: the smallest non-singleton
// color class, ties broken by (smallest) color value.  This rule is a
// pure function of the partition, hence isomorphism-invariant — which is
// what makes the set of leaves the search explores correspond under any
// graph isomorphism, and so the lexicographic-minimum leaf render a
// true canonical form.  Returns the chosen color, or -1 if the coloring
// is already discrete (never happens when called on a non-leaf).
int pick_cell(const std::vector<int>& color) {
  std::map<int, int> count;
  for (int const c : color)
    ++count[c];
  int chosen = -1;
  int best_size = 0;
  for (const auto& [c, s] : count) { // map iterates in ascending color order
    if (s < 2)
      continue;
    if (chosen < 0 || s < best_size) {
      chosen = c;
      best_size = s;
    }
  }
  return chosen;
}

// Individualization-refinement search state (plan §3.2 step 3).
// Inputs are held by pointer (not reference) so the struct stays a
// plain value type — clang-tidy gates reference data members.
struct SearchState {
  const ComplexGraph* g;
  const std::vector<std::vector<int>>* adj;
  const std::vector<int>* comp_vertex;
  int n_mol;
  long leaf_budget;      // remaining leaves before the §5 guard trips
  bool best_set = false; // false until the first leaf is rendered
  std::string best;      // lexicographically minimal leaf render so far
};

// Recurse: `color` is a WL-stable coloring.  At a leaf, render and keep
// the minimum.  Otherwise pick a target cell and, for each of its
// vertices, individualize that vertex (give it a fresh top color),
// re-refine, and recurse.  Individualization strictly splits the cell,
// so every path reaches a discrete (hence leaf) coloring in at most
// n_vert steps — termination needs no separate depth bound; the leaf
// budget guards only against a pathological branching factor.
void search(SearchState& st, const std::vector<int>& color) {
  if (is_leaf(*st.g, color, *st.comp_vertex, st.n_mol)) {
    std::string label = render_leaf(*st.g, color, *st.comp_vertex, st.n_mol);
    if (!st.best_set || label < st.best) {
      st.best = std::move(label);
      st.best_set = true;
    }
    --st.leaf_budget;
    return;
  }
  if (st.leaf_budget <= 0)
    return; // §5 safety valve — keep the best leaf already found

  int const target = pick_cell(color);
  int mx = 0;
  for (int const c : color)
    mx = std::max(mx, c);

  int const n_vert = static_cast<int>(color.size());
  for (int v = 0; v < n_vert; ++v) {
    if (color[v] != target)
      continue;
    if (st.leaf_budget <= 0)
      return;
    // Individualize v: a fresh color above every existing one.  refine
    // re-ranks, leaving v alone in the top class; the choice is applied
    // identically to isomorphic graphs, so the search trees correspond.
    std::vector<int> branch = color;
    branch[v] = mx + 1;
    refine(*st.adj, branch);
    search(st, branch);
  }
}

} // namespace

CanonForm canonicalize(const ComplexGraph& g) {
  int const n_mol = g.molecule_count();
  if (n_mol == 0)
    return {"", true};

  // --- Build the port graph ------------------------------------------------
  //
  // Vertices: [0, n_mol)        molecule vertices
  //           [n_mol, n_vert)   one per bonded component
  //
  // comp_vertex[gc] is the vertex id of bonded component `gc`, or -1
  // for an unbonded component (no vertex — it folds into its molecule's
  // color).
  std::vector<int> comp_vertex(g.components.size(), -1);
  std::vector<int> vertex_comp; // component-vertex index -> global comp index
  std::vector<int> comp_to_mol(g.components.size(), -1);
  for (int m = 0; m < n_mol; ++m) {
    const auto& mol = g.molecules[m];
    for (int i = 0; i < mol.n_comp; ++i) {
      int const gc = mol.first_comp + i;
      comp_to_mol[gc] = m;
      if (g.components[gc].partner >= 0) {
        comp_vertex[gc] = n_mol + static_cast<int>(vertex_comp.size());
        vertex_comp.push_back(gc);
      }
    }
  }
  int const n_vert = n_mol + static_cast<int>(vertex_comp.size());

  // Adjacency: molecule <-> bonded-component edges, then bond edges.
  std::vector<std::vector<int>> adj(n_vert);
  for (int cv = 0; cv < static_cast<int>(vertex_comp.size()); ++cv) {
    int const vid = n_mol + cv;
    int const mv = comp_to_mol[vertex_comp[cv]];
    adj[vid].push_back(mv);
    adj[mv].push_back(vid);
  }
  for (int cv = 0; cv < static_cast<int>(vertex_comp.size()); ++cv) {
    int const gc = vertex_comp[cv];
    int const partner = g.components[gc].partner;
    if (gc < partner) { // emit each bond edge once
      int const va = n_mol + cv;
      int const vb = comp_vertex[partner];
      adj[va].push_back(vb);
      adj[vb].push_back(va);
    }
  }

  // --- Initial colors ------------------------------------------------------
  //
  // Intern color strings through a sorted map, so the assigned integers
  // follow lexicographic string order — the integer colors are a
  // canonical, structure-derived ranking from round zero.
  std::vector<std::string> init_str(n_vert);
  for (int m = 0; m < n_mol; ++m)
    init_str[m] = molecule_color(g, g.molecules[m]);
  for (int cv = 0; cv < static_cast<int>(vertex_comp.size()); ++cv)
    init_str[n_mol + cv] = component_color(g.components[vertex_comp[cv]]);

  std::map<std::string, int> intern;
  for (const auto& s : init_str)
    intern.emplace(s, 0);
  {
    int next = 0;
    for (auto& [str, id] : intern)
      id = next++;
  }
  std::vector<int> color(n_vert);
  for (int v = 0; v < n_vert; ++v)
    color[v] = intern[init_str[v]];

  // --- 1-WL color refinement ----------------------------------------------
  refine(adj, color);

  // --- Fast path -----------------------------------------------------------
  //
  // If refinement alone discriminated the complex (plan §3.2 step 2),
  // the refined colors fix a unique isomorphism-invariant ordering;
  // render directly.  This is the overwhelmingly common case.
  if (is_leaf(g, color, comp_vertex, n_mol))
    return {render_leaf(g, color, comp_vertex, n_mol), /*fast_path=*/true};

  // --- Individualization-refinement (plan §3.2 step 3) ---------------------
  //
  // A genuinely symmetric complex (rings, homo-oligomers) survived
  // refinement.  Branch on a target cell, individualize each member,
  // re-refine, and recurse; the canonical label is the lexicographically
  // minimal leaf render.  This is a true canonical form — see search()
  // and pick_cell() for why the leaf set is isomorphism-invariant.
  SearchState st{&g, &adj, &comp_vertex, n_mol, kLeafBudget, false, std::string{}};
  search(st, color);
  return {st.best, /*fast_path=*/false};
}

std::string canonical_label(const ComplexGraph& g) { return canonicalize(g).label; }

} // namespace rulemonkey::canonical
