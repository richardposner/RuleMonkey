// Faithful port of NFsim's Sekar energy-rule expansion.
// See energy_expand.hpp for provenance and the governing algorithm.

#include "energy_expand.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <utility>

namespace rulemonkey {

namespace {

// Internal working record for a context condition extracted from the
// conditional energy patterns.  Mirrors NFcore::ContextCondition.
struct ContextCondition {
  std::string mol_type;
  int reactant_idx = 0;
  std::string comp_name;
  std::vector<int> gated_pattern_indices; // patterns that require this bond
};

} // namespace

std::vector<int> EnergyFunction::find_relevant_for_binding(const std::string& mol_type1,
                                                           const std::string& site1,
                                                           const std::string& mol_type2,
                                                           const std::string& site2) const {
  std::vector<int> relevant;
  for (int i = 0; i < static_cast<int>(patterns_.size()); ++i) {
    const EnergyPatternInfo& ep = patterns_[i];
    for (const auto& bond : ep.bonds) {
      const EpMolecule& m1 = ep.molecules[bond.mol1];
      const EpMolecule& m2 = ep.molecules[bond.mol2];
      const std::string& c1 = m1.components[bond.comp1].name;
      const std::string& c2 = m2.components[bond.comp2].name;

      bool const match_forward =
          (m1.type_name == mol_type1 && c1 == site1 && m2.type_name == mol_type2 && c2 == site2);
      bool const match_reverse =
          (m1.type_name == mol_type2 && c1 == site2 && m2.type_name == mol_type1 && c2 == site1);
      if (match_forward || match_reverse) {
        relevant.push_back(i);
        break;
      }
    }
  }
  return relevant;
}

std::vector<EnergyExpandedVariant> EnergyFunction::expand_binding(const std::string& mol_type1,
                                                                  const std::string& site1,
                                                                  const std::string& mol_type2,
                                                                  const std::string& site2) const {
  std::vector<EnergyExpandedVariant> variants;

  // Step 1: relevant patterns (contain the reaction-center bond).
  std::vector<int> const relevant = find_relevant_for_binding(mol_type1, site1, mol_type2, site2);

  if (relevant.empty()) {
    // No energy patterns overlap the reaction center → ΔG = 0 everywhere.
    EnergyExpandedVariant v;
    v.delta_g = 0.0;
    variants.push_back(std::move(v));
    return variants;
  }

  // Step 2: classify relevant patterns into "always" (only the reaction
  // center) and "conditional" (extra bound/stated context or a third
  // molecule type).
  std::vector<int> always_patterns;
  std::vector<int> conditional_patterns;
  for (int const pi : relevant) {
    const EnergyPatternInfo& ep = patterns_[pi];
    bool has_extra_context = false;
    for (const auto& mol : ep.molecules) {
      bool const is_reactant_type = (mol.type_name == mol_type1 || mol.type_name == mol_type2);
      if (!is_reactant_type) {
        has_extra_context = true; // a third molecule type → conditional
        break;
      }
      for (const auto& comp : mol.components) {
        if (mol.type_name == mol_type1 && comp.name == site1)
          continue; // reaction-center site
        if (mol.type_name == mol_type2 && comp.name == site2)
          continue;
        if (comp.is_bound || !comp.state_constraint.empty()) {
          has_extra_context = true;
          break;
        }
      }
      if (has_extra_context)
        break;
    }
    if (has_extra_context)
      conditional_patterns.push_back(pi);
    else
      always_patterns.push_back(pi);
  }

  // Base ΔG from always-matching patterns.
  double base_g = 0.0;
  for (int const pi : always_patterns)
    base_g += patterns_[pi].energy_value;

  // Build the symbolic ΔG (sum of the given patterns' source expressions, or
  // "0"), so the caller can synthesize a set_param-re-resolvable rate.
  auto build_dg_expr = [this](const auto& indices) -> std::string {
    std::string expr;
    for (int const pi : indices) {
      std::string term = patterns_[pi].energy_expr;
      if (term.empty())
        term = std::to_string(patterns_[pi].energy_value);
      if (!expr.empty())
        expr += "+";
      expr += "(" + term + ")";
    }
    return expr.empty() ? "0" : expr;
  };

  // Step 3: extract context conditions from conditional patterns.  A
  // condition is a bound, non-center component on a reactant-type molecule;
  // the constraint added to the template is "bound to anything" (matching
  // NFsim's addBoundComponent), so the specific partner is not recorded.
  std::map<std::pair<std::string, std::string>, ContextCondition> cond_map;
  std::vector<std::pair<std::string, std::string>> cond_order; // stable insertion order
  for (int const pi : conditional_patterns) {
    const EnergyPatternInfo& ep = patterns_[pi];
    for (const auto& mol : ep.molecules) {
      int reactant_idx = -1;
      if (mol.type_name == mol_type1)
        reactant_idx = 0;
      else if (mol.type_name == mol_type2)
        reactant_idx = 1;
      else
        continue;

      for (const auto& comp : mol.components) {
        if (reactant_idx == 0 && comp.name == site1)
          continue;
        if (reactant_idx == 1 && comp.name == site2)
          continue;
        if (!comp.is_bound)
          continue;

        auto key = std::make_pair(mol.type_name, comp.name);
        auto it = cond_map.find(key);
        if (it == cond_map.end()) {
          ContextCondition cc;
          cc.mol_type = mol.type_name;
          cc.reactant_idx = reactant_idx;
          cc.comp_name = comp.name;
          cond_map[key] = std::move(cc);
          cond_order.push_back(key);
        }
        cond_map[key].gated_pattern_indices.push_back(pi);
      }
    }
  }

  // Materialize conditions in stable insertion order (NFsim's std::map
  // iteration is by key; insertion order keeps our variant numbering
  // deterministic and independent of type-name collation — the resulting
  // rule *set* is identical either way).
  std::vector<ContextCondition> conditions;
  conditions.reserve(cond_order.size());
  for (const auto& key : cond_order) {
    ContextCondition& cc = cond_map[key];
    std::sort(cc.gated_pattern_indices.begin(), cc.gated_pattern_indices.end());
    cc.gated_pattern_indices.erase(
        std::unique(cc.gated_pattern_indices.begin(), cc.gated_pattern_indices.end()),
        cc.gated_pattern_indices.end());
    conditions.push_back(cc);
  }

  if (conditions.empty()) {
    // All relevant patterns are "always" — one variant, ΔG = base_g.
    EnergyExpandedVariant v;
    v.delta_g = base_g;
    v.delta_g_expr = build_dg_expr(always_patterns);
    variants.push_back(std::move(v));
    return variants;
  }

  // Step 4: enumerate all 2^n context combinations.
  int const n_cond = static_cast<int>(conditions.size());
  int const n_combinations = 1 << n_cond;
  for (int combo = 0; combo < n_combinations; ++combo) {
    std::set<int> active_patterns(always_patterns.begin(), always_patterns.end());

    EnergyExpandedVariant v;
    for (int ci = 0; ci < n_cond; ++ci) {
      bool const condition_met = ((combo >> ci) & 1) != 0;
      EnergyContextConstraint cc;
      cc.reactant_idx = conditions[ci].reactant_idx;
      cc.comp_name = conditions[ci].comp_name;
      cc.must_be_bound = condition_met;
      v.constraints.push_back(std::move(cc));
      if (condition_met)
        for (int const pi : conditions[ci].gated_pattern_indices)
          active_patterns.insert(pi);
    }

    double delta_g = 0.0;
    for (int const pi : active_patterns)
      delta_g += patterns_[pi].energy_value;

    v.delta_g = delta_g;
    v.delta_g_expr = build_dg_expr(active_patterns);
    variants.push_back(std::move(v));
  }

  return variants;
}

namespace {

std::string site_name(const std::string& type, const std::string& comp) {
  return type + "." + comp;
}

// Why expand_binding cannot price one energy pattern exactly at the reaction
// center mol_type1.site1-mol_type2.site2, or "" if it can (see
// EnergyFunction::inexact_reason).  Each loop only finds what is wrong; the
// message is built after it.
std::string pattern_inexact_reason(const EnergyPatternInfo& ep, const std::string& mol_type1,
                                   const std::string& site1, const std::string& mol_type2,
                                   const std::string& site2, const SiteGraph& sites) {
  std::string const what = "energy pattern " + (ep.text.empty() ? ep.id : ep.text);
  std::string const center_bond = site_name(mol_type1, site1) + "-" + site_name(mol_type2, site2);
  auto is_center_site = [&](const std::string& type, const std::string& comp) {
    return (type == mol_type1 && comp == site1) || (type == mol_type2 && comp == site2);
  };
  int const n_mol = static_cast<int>(ep.molecules.size());

  // Which components are bond ends, and whether the bonds join every
  // molecule into one graph.
  std::vector<std::vector<bool>> bonded(n_mol);
  for (int mi = 0; mi < n_mol; ++mi)
    bonded[mi].assign(ep.molecules[mi].components.size(), false);
  std::vector<int> root(n_mol);
  for (int mi = 0; mi < n_mol; ++mi)
    root[mi] = mi;
  auto find = [&root](int m) {
    while (root[m] != m)
      m = root[m] = root[root[m]];
    return m;
  };
  for (const auto& b : ep.bonds) {
    bonded[b.mol1][b.comp1] = true;
    bonded[b.mol2][b.comp2] = true;
    root[find(b.mol1)] = find(b.mol2);
  }
  bool connected = true;
  for (int mi = 1; mi < n_mol; ++mi)
    connected = connected && find(mi) == find(0);
  if (!connected)
    return what + " is not connected, so joining or parting any two complexes can change its count";

  // The center bond, and every other bond.
  int n_center = 0;
  int center_x = -1; // the center bond's molecules
  int center_y = -1;
  std::vector<EnergyPatternInfo::Bond> context;
  for (const auto& b : ep.bonds) {
    const EpMolecule& m1 = ep.molecules[b.mol1];
    const EpMolecule& m2 = ep.molecules[b.mol2];
    const std::string& c1 = m1.components[b.comp1].name;
    const std::string& c2 = m2.components[b.comp2].name;
    bool const fwd =
        m1.type_name == mol_type1 && c1 == site1 && m2.type_name == mol_type2 && c2 == site2;
    bool const rev =
        m1.type_name == mol_type2 && c1 == site2 && m2.type_name == mol_type1 && c2 == site1;
    if (fwd || rev) {
      ++n_center;
      center_x = b.mol1;
      center_y = b.mol2;
    } else {
      context.push_back(b);
    }
  }

  // The first component asking for something the expansion cannot read:
  // with no center bond, a center site required free or bound to anything
  // (which a connected pattern needs for its count to change at all); with
  // one, any state and any site required free.  The `!+` components are
  // counted, and read below.
  enum class Ask { None, State, Free, Constrains };
  Ask ask = Ask::None;
  const EpMolecule* ask_mol = nullptr;
  const EpComponent* ask_comp = nullptr;
  int n_unnamed = 0;
  const EpMolecule* unnamed_mol = nullptr;
  const EpComponent* unnamed_comp = nullptr;
  for (int mi = 0; mi < n_mol; ++mi) {
    const EpMolecule& mol = ep.molecules[mi];
    for (std::size_t ci = 0; ci < mol.components.size(); ++ci) {
      const EpComponent& c = mol.components[ci];
      bool const unnamed = c.is_bound && !bonded[mi][ci];
      if (n_center != 0 && unnamed) {
        ++n_unnamed;
        unnamed_mol = &mol;
        unnamed_comp = &c;
      }
      if (ask != Ask::None)
        continue;
      if (n_center == 0) {
        if (is_center_site(mol.type_name, c.name) && (c.is_free || unnamed))
          ask = Ask::Constrains;
      } else if (!c.state_constraint.empty()) {
        ask = Ask::State;
      } else if (c.is_free) {
        ask = Ask::Free;
      }
      if (ask != Ask::None) {
        ask_mol = &mol;
        ask_comp = &c;
      }
    }
  }
  if (n_center == 0) {
    if (ask == Ask::None)
      return ""; // the bond cannot change its count
    return what + " constrains " + site_name(ask_mol->type_name, ask_comp->name) +
           " without naming the " + center_bond + " bond";
  }
  if (ask == Ask::State)
    return what + " requires " + site_name(ask_mol->type_name, ask_comp->name) + " in state " +
           ask_comp->state_constraint;
  if (ask == Ask::Free)
    return what + " requires " + site_name(ask_mol->type_name, ask_comp->name) + " free";

  // A component listed as `!?` asks nothing, but BNG2 counts its matches:
  // k of them, interchangeable, on the m same-named components the rest of
  // the pattern leaves give m!/(m-k)! matches over k! automorphisms, so
  // C(m, k) per match of the rest.  The expansion counts one.
  const EpMolecule* loose_mol = nullptr;
  std::string loose_name;
  long loose_ways = 1;
  for (int mi = 0; mi < n_mol && loose_ways == 1; ++mi) {
    const EpMolecule& mol = ep.molecules[mi];
    std::map<std::string, std::pair<int, int>> uses; // name -> (`!?`, other)
    for (const auto& c : mol.components) {
      bool const asks_nothing = !c.is_bound && !c.is_free && c.state_constraint.empty();
      ++(asks_nothing ? uses[c.name].first : uses[c.name].second);
    }
    for (const auto& [name, use] : uses) {
      if (use.first == 0)
        continue;
      auto it = sites.repeated.find({mol.type_name, name});
      int const declared = it == sites.repeated.end() ? 1 : it->second;
      int const m = declared - use.second;
      long ways = 1;
      for (int k = 0; k < use.first; ++k)
        ways = ways * std::max(m - k, 0) / (k + 1);
      if (ways != 1) {
        loose_mol = &mol;
        loose_name = name;
        loose_ways = ways;
        break;
      }
    }
  }
  if (loose_ways != 1)
    return what + " can place its " + site_name(loose_mol->type_name, loose_name) + " (!?) " +
           std::to_string(loose_ways) + " ways, which the expansion counts as one";

  if (n_center > 1)
    return what + " names the " + center_bond + " bond more than once";
  if (context.size() > 1)
    return what + " names more than one bond besides the " + center_bond + " bond";
  if (context.size() + n_unnamed > 1)
    return what + " requires more than one site besides the " + center_bond + " bond bound";

  if (n_unnamed == 1) {
    // A lone `!+`, on one of the two center molecules, is the expansion's
    // own condition, "this site is bound to anything", when the site is not
    // named like the reaction site and its molecule type declares one.
    SiteKey const key{unnamed_mol->type_name, unnamed_comp->name};
    std::string const where = site_name(key.first, key.second);
    if (is_center_site(key.first, key.second))
      return what + " names a second " + where + " besides the reaction site";
    if (sites.repeated.count(key) != 0)
      return what + " requires " + where + " bound, but " + key.first + " declares more than one " +
             key.second;
    return "";
  }
  if (context.empty())
    return ""; // the center bond alone

  // One more bond, which must run from a center molecule to a molecule of a
  // third type.
  const auto& b = context.front();
  bool const in1 = b.mol1 == center_x || b.mol1 == center_y;
  bool const in2 = b.mol2 == center_x || b.mol2 == center_y;
  if (in1 && in2)
    return what + " closes a ring with the " + center_bond + " bond";
  const EpMolecule& near_mol = ep.molecules[in1 ? b.mol1 : b.mol2];
  const EpMolecule& far_mol = ep.molecules[in1 ? b.mol2 : b.mol1];
  SiteKey const near{near_mol.type_name, near_mol.components[in1 ? b.comp1 : b.comp2].name};
  SiteKey const far{far_mol.type_name, far_mol.components[in1 ? b.comp2 : b.comp1].name};
  std::string const near_site = site_name(near.first, near.second);
  if (far.first == mol_type1 || far.first == mol_type2)
    return what + " names a second " + far.first + ", a molecule type of the reaction center";
  if (is_center_site(near.first, near.second))
    return what + " names a second " + near_site + " besides the reaction site";
  if (sites.repeated.count(near) != 0)
    return what + " requires " + near_site + " bound, but " + near.first +
           " declares more than one " + near.second;
  auto it = sites.partners.find(near);
  if (it == sites.partners.end())
    return "";
  auto other = std::find_if(it->second.begin(), it->second.end(),
                            [&far](const SiteKey& p) { return p != far; });
  if (other == it->second.end())
    return "";
  // The pattern's own pairing may never occur at all; "also" only if it can.
  std::string const binds = it->second.count(far) != 0 ? " can also bind " : " binds ";
  return what + " requires " + near_site + " bound to " + site_name(far.first, far.second) +
         ", but " + near_site + binds + site_name(other->first, other->second);
}

} // namespace

std::string EnergyFunction::inexact_reason(const std::string& mol_type1, const std::string& site1,
                                           const std::string& mol_type2, const std::string& site2,
                                           const SiteGraph& sites) const {
  for (const EnergyPatternInfo& ep : patterns_) {
    std::string why = pattern_inexact_reason(ep, mol_type1, site1, mol_type2, site2, sites);
    if (!why.empty())
      return why;
  }
  return "";
}
} // namespace rulemonkey
