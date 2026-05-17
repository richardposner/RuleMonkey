#pragma once

// Canonical complex labeling — the pure species-identity primitive.
//
// Design & decisions: dev/canonical_labeling_plan_2026_05_16.md.
// Issue #9 §2 (`.species` output) is the first consumer; partial
// scaling (plan §7.2) is the consumer that shapes the design.
//
// This header is the PURE CORE only (plan decision #3): a value type
// `ComplexGraph` and the pure function `canonical_label`.  It has NO
// dependency on AgentPool / Engine internals, so it is callable from
// both the on-demand batch sweep and the (future) cached-incremental
// layer (decision #6), and is unit-testable with hand-built graphs and
// no simulator.
//
// `extract_complex` (engine -> ComplexGraph) is deliberately NOT
// declared here: AgentPool is TU-local to engine.cpp by design, so the
// engine->graph bridge lives with the engine.  See the plan handoff /
// the step-0 report.
//
// SCOPE as of this module's first landing (plan §6 steps 1-2):
// `canonical_label` implements 1-WL color refinement + the
// all-distinct-colors fast path.  Individualization-refinement for
// genuinely symmetric complexes (rings, homo-oligomers — plan §3.2
// step 3) is NOT yet implemented; for such inputs the returned label
// is deterministic per input but not guaranteed isomorphism-invariant.
// `canonicalize().fully_refined` reports which path was taken.

#include <string>
#include <utility>
#include <vector>

namespace rulemonkey::canonical {

// ---------------------------------------------------------------------------
// ComplexGraph — the vertex-colored port graph of one connected complex.
//
// Value type only: molecule / component / bond structure, no engine
// reference.  Hand-buildable for tests; produced by the engine-side
// `extract_complex` for live pools.
//
// Components are stored in one flat array, grouped by molecule and in
// each molecule's BNGL declaration order (the order the renderer
// emits).  A bond is recorded symmetrically as a `partner` global
// component index on both endpoints (-1 = free), mirroring the
// engine's ComponentInstance representation.
// ---------------------------------------------------------------------------
struct ComplexGraph {
  struct Component {
    std::string name;
    std::string state; // "" = stateless component (no `~state` rendered)
    int partner = -1;  // global component index of bonded partner; -1 = free
  };
  struct Molecule {
    std::string type_name;
    int first_comp = 0; // index of this molecule's first component
    int n_comp = 0;     // component count
  };

  std::vector<Molecule> molecules;
  std::vector<Component> components; // flat, grouped by molecule

  // --- Builder -------------------------------------------------------------

  // Append a molecule with its components in declaration order.
  // `comps` is a list of (name, state) pairs; pass "" for a stateless
  // component.  Returns the new molecule's index.
  int add_molecule(const std::string& type_name,
                   const std::vector<std::pair<std::string, std::string>>& comps);

  // Bond two components, each addressed as (molecule index, local
  // component index within that molecule).  A self-bond (same molecule,
  // two different components) is permitted.
  void add_bond(int mol_a, int local_a, int mol_b, int local_b);

  // Global component index of molecule `mol`'s local component `local`.
  int global_comp(int mol, int local) const { return molecules[mol].first_comp + local; }
  int molecule_count() const { return static_cast<int>(molecules.size()); }
  int component_count() const { return static_cast<int>(components.size()); }
};

// ---------------------------------------------------------------------------
// Canonicalization
// ---------------------------------------------------------------------------

// Result of canonicalizing one complex.
struct CanonForm {
  std::string label;  // canonical BNGL species string (also the dedup key)
  bool fully_refined; // true  -> 1-WL fully discriminated every molecule;
                      //          `label` is a true canonical form.
                      // false -> a symmetric color class survived
                      //          refinement; `label` used an incomplete
                      //          tie-break and is NOT guaranteed
                      //          isomorphism-invariant until plan §3.2
                      //          step 3 (individualization) lands.
};

// Canonicalize a complex: 1-WL color refinement, then render the
// canonical BNGL string by walking molecules in refined-color order.
// Pure function of the graph — see `CanonForm::fully_refined` for the
// current correctness boundary.
CanonForm canonicalize(const ComplexGraph& g);

// Convenience wrapper: the canonical BNGL string only.  This is the
// signature pinned by the plan (§4); both calling modes use it.
std::string canonical_label(const ComplexGraph& g);

} // namespace rulemonkey::canonical
