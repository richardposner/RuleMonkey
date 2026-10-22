// Exact Arrhenius free-energy changes: the parts of energy_match.hpp that
// are not templates.  See that header for what is counted and why.

#include "energy_match.hpp"

#include <numeric>

namespace rulemonkey {

namespace {

// A pattern read as a graph, so that matching it on itself counts its
// automorphisms the way BioNetGen does.  Components are numbered flat, in
// pattern order.  A `!+` is bound to something off the graph, and a `!?` to
// something unknown, which matches only a `!?`.
class PatternView {
public:
  explicit PatternView(const EmPattern& p) : p_(p) {
    for (int m = 0; m < static_cast<int>(p.mols.size()); ++m) {
      first_.push_back(static_cast<int>(mol_of_.size()));
      for (std::size_t c = 0; c < p.mols[m].comps.size(); ++c)
        mol_of_.push_back(m);
    }
  }
  int type_of(int mol) const { return p_.mols[mol].type; }
  int mol_of(int comp) const { return mol_of_[comp]; }
  int name_of(int comp) const { return at(comp).name; }
  int state_of(int comp) const { return at(comp).state; }
  int partner_of(int comp) const {
    const EmComponent& c = at(comp);
    switch (c.bond) {
    case EmComponent::Bond::To:
      return first_[c.to_mol] + c.to_comp;
    case EmComponent::Bond::Bound:
      return kEmPhantom;
    case EmComponent::Bond::Free:
      return -1;
    case EmComponent::Bond::Any:
      return kEmUnknown;
    }
    return kEmUnknown;
  }
  void comps_named(int mol, int name, std::vector<int>& out) const {
    const auto& comps = p_.mols[mol].comps;
    for (int c = 0; c < static_cast<int>(comps.size()); ++c)
      if (comps[c].name == name)
        out.push_back(first_[mol] + c);
  }
  int first(int mol) const { return first_[mol]; }

private:
  const EmComponent& at(int comp) const {
    int const m = mol_of_[comp];
    return p_.mols[m].comps[comp - first_[m]];
  }
  const EmPattern&
      p_; // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members): lives inside one count
  std::vector<int> first_;
  std::vector<int> mol_of_;
};

// The molecules reached from `start` by the pattern's bonds, leaving out
// the bond (cut_mol, cut_comp) in both directions.
std::vector<int> reach(const EmPattern& p, int start, int cut_mol, int cut_comp) {
  const EmComponent& cut = p.mols[cut_mol].comps[cut_comp];
  std::vector<char> seen(p.mols.size(), 0);
  std::vector<int> order{start};
  seen[start] = 1;
  for (std::size_t k = 0; k < order.size(); ++k) {
    int const m = order[k];
    const auto& comps = p.mols[m].comps;
    for (int ci = 0; ci < static_cast<int>(comps.size()); ++ci) {
      const EmComponent& ec = comps[ci];
      if (ec.bond != EmComponent::Bond::To)
        continue;
      bool const is_cut =
          (m == cut_mol && ci == cut_comp) || (m == cut.to_mol && ci == cut.to_comp);
      if (is_cut || seen[ec.to_mol])
        continue;
      seen[ec.to_mol] = 1;
      order.push_back(ec.to_mol);
    }
  }
  return order;
}

// The sub-pattern on `mols` (indices into p), with the cut component's bond
// dropped; `root` is translated to its index in the sub-pattern.
EmPattern sub_pattern(const EmPattern& p, const std::vector<int>& mols, int cut_mol, int cut_comp,
                      int& root_mol) {
  std::vector<int> index(p.mols.size(), -1);
  for (int k = 0; k < static_cast<int>(mols.size()); ++k)
    index[mols[k]] = k;
  EmPattern side;
  side.text = p.text;
  for (int const m : mols) {
    EmMolecule mol = p.mols[m];
    for (int ci = 0; ci < static_cast<int>(mol.comps.size()); ++ci) {
      EmComponent& ec = mol.comps[ci];
      if (m == cut_mol && ci == cut_comp) {
        ec.bond = EmComponent::Bond::Any;
        ec.to_mol = -1;
        ec.to_comp = -1;
      } else if (ec.bond == EmComponent::Bond::To) {
        ec.to_mol = index[ec.to_mol];
      }
    }
    side.mols.push_back(std::move(mol));
  }
  root_mol = index[cut_mol];
  return side;
}

double falling(int n, int k) {
  double out = 1.0;
  for (int i = 0; i < k; ++i)
    out *= std::max(n - i, 0);
  return out;
}

} // namespace

bool em_connected(const EmPattern& p) {
  if (p.mols.empty())
    return true;
  std::vector<int> root(p.mols.size());
  std::iota(root.begin(), root.end(), 0);
  auto find = [&root](int m) {
    while (root[m] != m)
      m = root[m] = root[root[m]];
    return m;
  };
  for (const auto& mol : p.mols)
    for (const auto& ec : mol.comps)
      if (ec.bond == EmComponent::Bond::To)
        root[find(static_cast<int>(&mol - p.mols.data()))] = find(ec.to_mol);
  for (std::size_t m = 1; m < p.mols.size(); ++m)
    if (find(static_cast<int>(m)) != find(0))
      return false;
  return true;
}

double em_automorphisms(const EmPattern& p) {
  if (p.mols.empty() || p.mols[0].comps.empty())
    return 1.0; // one molecule, nothing to permute
  PatternView const view(p);
  // Every self-map puts molecule 0's first component somewhere: sum over
  // where.
  long total = 0;
  std::vector<int> targets;
  for (int m = 0; m < static_cast<int>(p.mols.size()); ++m)
    if (p.mols[m].type == p.mols[0].type)
      view.comps_named(m, p.mols[0].comps[0].name, targets);
  for (int const t : targets)
    total += em_count_rooted(p, view, 0, 0, t);
  return total > 0 ? static_cast<double>(total) : 1.0;
}

bool em_cut(const EmPattern& p, int pm, int pc, EmPattern& side0, int& root_mol0, int& root_comp0,
            EmPattern& side1, int& root_mol1, int& root_comp1) {
  const EmComponent& bond = p.mols[pm].comps[pc];
  int const pm2 = bond.to_mol;
  int const pc2 = bond.to_comp;
  std::vector<int> const mols0 = reach(p, pm, pm, pc);
  for (int const m : mols0)
    if (m == pm2)
      return false; // the bond is on a cycle
  std::vector<int> const mols1 = reach(p, pm2, pm, pc);
  side0 = sub_pattern(p, mols0, pm, pc, root_mol0);
  root_comp0 = pc;
  side1 = sub_pattern(p, mols1, pm2, pc2, root_mol1);
  root_comp1 = pc2;
  return true;
}

double em_rooted_bound(const EmPattern& side, int pm, int pc, const ArrheniusModel& am) {
  // Reach the molecules from the root as the matcher does, recording the
  // component each is reached through: the root's own, and for the rest the
  // near end of the bond that reaches it.  That component's image is fixed;
  // the molecule's other listed components choose among the type's
  // same-named slots left, k of them on n slots in at most n!/(n-k)! ways.
  std::vector<int> fixed_comp(side.mols.size(), -1);
  std::vector<int> order{pm};
  fixed_comp[pm] = pc;
  for (std::size_t k = 0; k < order.size(); ++k) {
    const auto& comps = side.mols[order[k]].comps;
    for (const auto& ec : comps)
      if (ec.bond == EmComponent::Bond::To && fixed_comp[ec.to_mol] < 0) {
        fixed_comp[ec.to_mol] = ec.to_comp;
        order.push_back(ec.to_mol);
      }
  }
  double bound = 1.0;
  for (int m = 0; m < static_cast<int>(side.mols.size()); ++m) {
    const EmMolecule& mol = side.mols[m];
    std::map<int, int> listed; // name -> components listed
    for (const auto& ec : mol.comps)
      ++listed[ec.name];
    const auto& slots = am.slot_names[mol.type];
    for (const auto& [name, k] : listed) {
      int const n = static_cast<int>(std::count(slots.begin(), slots.end(), name));
      int const fixed = (fixed_comp[m] >= 0 && mol.comps[fixed_comp[m]].name == name) ? 1 : 0;
      bound *= falling(n - fixed, k - fixed);
    }
  }
  return bound;
}

} // namespace rulemonkey
