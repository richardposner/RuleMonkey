/*
** Copyright 2008 Northern Arizona University and
** The Translational Genomics Research Institute
**
**  This file is part of Rulemonkey.
**
**  Rulemonkey is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  Rulemonkey is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with Rulemonkey.  If not, see <http://www.gnu.org/licenses/>.
**
*/


#include <assert.h>

#include "dynrules.h"
#include "dyncolors.h"
#include "dynmaps.h"
#include "rule.h"
#include "constants.h"
#include "naututil.h"
#include "output.h"
#include "perm_set.h"
#include "world.h"
extern world_t *world;

void suppress_warnings9() {SUPPRESS_WARNINGS;}

/* Local function prototypes. */
static int rule_output_raw(const output_type ptype, const int level,
                           rule_t *rule, dyncolors_t *dyncolors,
                           const char *prefix, boolean *used);
static void rule_ilist_output(const output_type ptype, const int level,
                                   const char *label,
                                   ilist_t *ilist);
static void rule_permute_cleanup(rule_t *molecule_rule);
static int rule_permute_all(specie_t *specie,
                            dynrule_aggregates_t *dynrule_aggregates,
                            rule_t *main_rule,
                            rule_t *rules_to_permute, int *pcount);
static int rule_permute_check(specie_t *specie,
                              dynrule_aggregates_t *dynrule_aggregates,
                              rule_t *rule, int *pcount);
static int rule_molecule_permute_check(specie_t *specie, rule_t *molocule_rule, 
                                       int *pcount);
static int rule_component_permute_check(specie_t *specie,
                                        rule_t *component_rule,
                                        int *pcount);
static int rule_molecule_check(specie_t *specie, rule_t *molecule_rule,
                               int *pcount, const int update);
static int rule_component_check(specie_t *specie, rule_t *component_rule,
                                const int molecule_node,
                                int *pcount, const int update);
static int rule_component_connection_check
             (specie_t *specie, const int molecule_node,
              rule_t *component_connection_rule,
              set *component_connection_row, int *pcount,
              const int update);
static int rule_attribute_check(int node, specie_t *specie, rule_t *rule);
static int rule_unique_molecule_node_check(rule_t *main_rule);
int rule_check_raw(specie_t *specie, dynrule_aggregates_t *dynrule_aggregates,
                   rule_t *rule, int *pcount, dynrules_t *dynrules,
                   dyncolors_t *dyncolors);
static boolean rule_count_check(rule_t *rule, const int count);
static int rule_perm_set_add_molecule(const int perm_base,
                                      rule_t *molecule_rule,
                                      const int molecule_node);
static int rule_perm_set_create(rule_t *molecule_rule);

/**
 * Handles or rules.
 */
int
rule_check(specie_t *specie, dynrule_aggregates_t *dynrule_aggregates,
           rule_t *rule, int *pcount, dynrules_t *dynrules,
           dyncolors_t *dyncolors)
{
  int count;

  *pcount = 0;

  while (NULL != rule)
  {
    int status;

    status = rule_check_raw(specie, dynrule_aggregates, rule, &count, dynrules,
                            dyncolors);
    if (FALSE == status) return FALSE;

    // Calculate sum of all or rules, required for observations.
    *pcount += count;

    rule = rule->or_rule;
  }

  return TRUE;
}

/**
 * Find all possible rule matches for the given specie.
 * First, each component rule is assigned a list of nodes that match
 * all the preconditions.
 * Second, all possible permutations of component rule to node mappings are
 * tested, and the valid ones are added to the rule aggregate.
 */
int
rule_check_raw(specie_t *specie, dynrule_aggregates_t *dynrule_aggregates,
               rule_t *rule, int *pcount, dynrules_t *dynrules,
               dyncolors_t *dyncolors)
{
  int count;
  int molecule_count;
  rule_t *molecule_rule;
  int status;

  molecule_count = 0;

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "# Checking rule %d, specie %d.\n",
                 rule->index, specie->id + 1);
#endif

  // Loop through molecule_rules to perform initial matching.
  for(molecule_rule = rule;
      NULL != molecule_rule;
      molecule_rule = molecule_rule->and_rule)
  {
#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "Checking molecule rule %d, specie %d.\n",
                 molecule_rule->index, specie->id + 1);
#endif

    // Check if molecule satisfies rule.
    status = rule_molecule_check(specie, molecule_rule, &count, FALSE);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }

    if (FALSE == rule_count_check(molecule_rule, count))
    {
      // All molecules must match for rule to be successful.
      *pcount = 0;

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "Rule %d, specie %d (molecule_rule %d) failed.\n",
                 rule->index, specie->id + 1, molecule_rule->index);
#endif

      // No error occured.
      return TRUE;
    }
  }

  // Successfully matched, so collect all component rule to node mappings.
  for(molecule_rule = rule;
      NULL != molecule_rule;
      molecule_rule = molecule_rule->and_rule)
  {
    // Check if molecule satisfies rule.
    status = rule_molecule_check(specie, molecule_rule, &count, TRUE);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }

    if (FALSE == rule_count_check(molecule_rule, count))
    {
      // Node mapping failed, so rules didn't match after all.

      // Just break out of loop so that rule_permute_cleanup gets called below.
      break;
    }
#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "molecule_check: Rule %d, molecule_rule %d, "
                 "specie %d matched %d times.\n",
                 rule->index, molecule_rule->index, specie->id + 1, count);
#endif
  }

  if (count != 0)
  {
    // Step through each component permutation, checking
    // component-component connections and adding to rule aggregates
    // when appropriate.
    count = 0;
#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "About to call rule_permute_all (rule:%d, specie:%d).\n",
                 rule->index, specie->id + 1);
    DYNRULES_OUTPUT_RULE(OUTPUT_DEBUG, 5, rule, dynrules, dyncolors);
#endif
    status = rule_permute_all(specie, dynrule_aggregates,
                              rule, rule, &count);
    if (FALSE == status)
    {
      return FALSE;
    }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "permute_all: Rule %d, specie %d matched %d times.\n",
                 rule->index, specie->id + 1, count);
#endif
  }
  else
  {
#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "Rule %d, specie %d failed while saving node maps.\n",
                 rule->index, specie->id + 1);
#endif
  }

  // Remove node maps and permute objects to reset rules.
  rule_permute_cleanup(rule);

  // Return number of permutations found.
  *pcount = count;

  return TRUE;
}

/**
 * Reset node maps on all component rules.  Node maps are only used in
 * component rules, so don't worry about any others.
 * Reset component permutations on all rules.  Component rules simply
 * have a pointer to the real perm_set, which are assigned to the
 * molecule rules.
 */
static void
rule_permute_cleanup(rule_t *rule)
{
  rule_t *molecule_rule;

  for (molecule_rule = rule;
       molecule_rule != NULL;
       molecule_rule = molecule_rule->and_rule)
  {
    rule_t *component_rule;

    for (component_rule = molecule_rule->adjacent_rule;
         component_rule != NULL;
         component_rule = component_rule->and_rule)
    {
      // Reset node map on component rule.
      component_rule->node_map->count = 0;

      // Remove pointer to molecule perm_set.
      component_rule->perm_set = NULL;
      component_rule->perm_set_index = 0;
    }

    // Destroy component permutation object.
    if (NULL != molecule_rule->perm_set)
    {
      perm_set_destroy(molecule_rule->perm_set);
      molecule_rule->perm_set = NULL;
      molecule_rule->perm_set_index = 0;
    }
  }
}


/**
 * Select each possible component connection in turn, and attempt to
 * complete the rule check on the given specie.
 */
static int
rule_permute_all(specie_t *specie, dynrule_aggregates_t *dynrule_aggregates,
                 rule_t *main_rule, rule_t *rules_to_permute, int *pcount)
{
  perm_set_t *perm_set;
  int set_count;
  int set_size;
  int i;
  int status;

  perm_set = rules_to_permute->perm_set;
  set_size = perm_set->set_size;
  set_count = perm_set_count(perm_set);

  // Loop through each permutation of component rules to node mappings from
  // current molecule.
  for (i = 0; i < set_count; i++)
  {
    perm_set->current_perm = perm_set_get_existing(i, perm_set);

    if (0 == perm_set_unique_check(set_size, perm_set->current_perm))
    {
      // Single node assigned to more than one component rule,
      // which is not physically possible, so ignore this extra permutation.
      continue;
    }

    if (NULL != rules_to_permute->and_rule)
    {
      // Continue applying permutations.
      status = rule_permute_all(specie, dynrule_aggregates, main_rule,
                                rules_to_permute->and_rule, pcount);
      if (FALSE == status)
      {
        return FALSE;
      }
    }
    else
    {
      // Everything narrowed down to one permutation.

      // Check that all molecule rules have unique nodes assigned.
      status = rule_unique_molecule_node_check(main_rule);
      if (FALSE == status)
      {
        // Single node assigned to more than one molecule rule,
        // which is not physically possible, so ignore this extra permutation.
#ifdef RULE_CHECK_DEBUG
        DEBUG_PRINTF(5, "Main rule %d (specie %d) failed "
                     "unique molecule node check.\n",
                     main_rule->index, specie->id + 1);
#endif

        continue;
      }
      
      // Continue with final rule check.
      status = rule_permute_check(specie, dynrule_aggregates,
                                  main_rule, pcount);
      if (FALSE == status)
      {
        return FALSE;
      }
    }
  }

  return TRUE;
}

/**
 * Go through all molecule rules and verify that each has a unique node
 * assigned to it.
 */
static int
rule_unique_molecule_node_check(rule_t *main_rule)
{
  for (;
       main_rule != NULL;
       main_rule = main_rule->and_rule)
  {
    rule_t *current_rule;
    int main_node;

    assert(0 == main_rule->perm_set_index);
    main_node = main_rule->perm_set->current_perm[main_rule->perm_set_index];
    for (current_rule = main_rule->and_rule;
         current_rule != NULL;
         current_rule = current_rule->and_rule)
    {
      if (main_node ==
          current_rule->perm_set->current_perm[current_rule->perm_set_index])
      {
        // single node assigned to more than one component.
        return FALSE;
      }
    }
  }

  return TRUE;
}

/**
 * Check if all current component connections match rule requirements,
 * adding an entry into rule aggregate if it is successful.
 * Start by looping through each molecule rule.
 */
static int
rule_permute_check(specie_t *specie, dynrule_aggregates_t *dynrule_aggregates,
                   rule_t *rule, int *pcount)
{
  int status;
  int count;
  rule_t *molecule_rule;

  count = 0;

  // Loop through molecule_rules to verify all component-component
  // connection permutations.
  for (molecule_rule = rule;
       NULL != molecule_rule;
       molecule_rule = molecule_rule->and_rule)
  {
    // Check if molecule satisfies rule.
    status = rule_molecule_permute_check(specie, molecule_rule,
                                          &count);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(8, "molecule_permute_check: molecule_rule %d, rule %d, "
                 "specie %d, count %d\n",
                 molecule_rule->index, rule->index, specie->id + 1, count);
#endif
    if (count <= 0)
    {
      // A rule failed to match.
      // No errors occured.
      return TRUE;
    }
  }

  // Don't try to add to dynrule_aggregates if it isn't present.
  if (NULL != dynrule_aggregates)
  {
    // Current permutation is completely verified, so add all component rules
    // from all molecule rules to rule aggregate.

    // First allocate new aggregate.
    status = dynrule_aggregates_rule_aggregate_create(dynrule_aggregates);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }

    // Loop through all molecule rules to apply all component rules.
    for (molecule_rule = rule;
         NULL != molecule_rule;
         molecule_rule = molecule_rule->and_rule)
    {
      rule_t *component_rule;

      for (component_rule = molecule_rule->adjacent_rule;
           NULL != component_rule;
           component_rule = component_rule->and_rule)
      {
        int matched_index;
        int matched_node;

        matched_index = component_rule->perm_set_index;
        matched_node = component_rule->perm_set->current_perm[matched_index];
        status = dynrule_aggregates_add_rule(dynrule_aggregates,
                                    component_rule, matched_node);
        if (FALSE == status)
        {
          // Problem while applying rule.
          assert(0);
          return FALSE;
        }
      }
    }

    // Now that we have a complete and valid aggregate, make sure it doesn't
    // duplicate a previously added aggregate.
    status = dynrule_aggregates_check_current_unique(dynrule_aggregates);
    if (1 == status)
    {
      // Aggregate was unique, now check that product count is correct.
#ifdef RULE_CHECK_DEBUG
      DEBUG_PRINTF(8, "About to call "
                   "dynrule_aggregates_check_product_count.\n");
      DEBUG_PRINTF(8, "dynrule aggregates:.\n");
      DYNRULE_AGGREGATES_OUTPUT(OUTPUT_DEBUG, 8, dynrule_aggregates);
#endif
      status = dynrule_aggregates_check_product_count
                 (dynrule_aggregates, specie);
    }

    if (0 == status)
    {
      // Either current aggregate is not unique or product count
      // is not correct.
      // New rule aggregate has already beeen removed.
      // Incomplete match, no errors occured.
      return TRUE;
    }
  }

  // Current permutation matched rules completely.
  // Increment number of permutations found.
  // If component rules are present, count will always be 1, so each
  // valid permutation will only be counted once.
  // If component rules are not present, count will be set to the number
  // of designated molecules in current specie.
  // This means that for a valid result, either specify a single molecule
  // with no components, or give components for all molecules.
  *pcount += count;

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(4, "Matched: rule %d, specie %d: ", rule->index, specie->id + 1);
  for (molecule_rule = rule;
       NULL != molecule_rule;
       molecule_rule = molecule_rule->and_rule)
  {
    int set_size;
    int i;
    int *pnode;

    set_size = molecule_rule->perm_set->set_size;
    // Don't output set if no component permutations present.
    // (if set_size is 1, only the last molecule examined will be present.
    if (set_size > 1)
    {
      DEBUG_PRINTF(4, "[");

      // Output current permutation.
      for (i = 0, pnode = molecule_rule->perm_set->current_perm;
           i < set_size;
           i++, pnode++)
      {
        if (i != 0)
        {
          DEBUG_PRINTF(4, ",");
        }
        DEBUG_PRINTF(4, "%d", *pnode);
      }
      DEBUG_PRINTF(4, "]");
    }
  }
  DEBUG_PRINTF(4, "\n");
#endif

  return TRUE;
}

/**
 * Loop through each component rule of the given molecule rule, performing
 * checks on specific permutation.
 */
static int
rule_molecule_permute_check(specie_t *specie, rule_t *molecule_rule, 
                            int *pcount)
{
  int count;
  rule_t *component_rule;

  *pcount = 0;

  if (NULL == molecule_rule->adjacent_rule)
  {
    // No components to match, so simply find out how many molecules are in
    // current specie, as they all match.
    // Need to provide sum now because no components are present to permute.

    // Count all molecules in current specie.
    *pcount = 1;

    return TRUE;
  }

  // Loop through component_rules.
  for (component_rule = molecule_rule->adjacent_rule;
       component_rule != NULL;
       component_rule = component_rule->and_rule)
  {
    int status;

    status = rule_component_permute_check(specie, component_rule,
                                          &count);
    if (FALSE == status) return FALSE;

    if (0 == count)
    {
      // Component check didn't pass.
      *pcount = 0;

      // Error did not occur.
      return TRUE;
    }
  }

  // All components checked out, so full molecule rule passes.
  *pcount = 1;

  return TRUE;
}

/**
 * For the given component rule, verify each and every component-component
 * connection rule.
 */
static int
rule_component_permute_check(specie_t *specie, rule_t *component_rule,
                             int *pcount)
{
  int component_node;
  set *component_row;
  int component_index;
  rule_t *adjacent_rule;
  int adjacent_index;
  int adjacent_node;

  *pcount = 0;

  if (NULL == component_rule->adjacent_rule)
  {
    // No component-component connections.  We've already verified these
    // are valid, so no need to check again.
    (*pcount)++;

    return TRUE;
  }

  component_index = component_rule->perm_set_index;
  component_node = component_rule->perm_set->current_perm[component_index];
  component_row = GRAPHROW(specie->cg->g.a, component_node, specie->cg->m);
  adjacent_rule = component_rule->adjacent_rule;

  adjacent_index = adjacent_rule->perm_set_index;
  // No further checking needed if wildcard rule.
  if (COLOR_WILDCARD != adjacent_rule->color_value)
  {
    adjacent_node = adjacent_rule->perm_set->current_perm[adjacent_index];
    if (!ISELEMENT(component_row, adjacent_node))
    {
      int curlen;

      curlen = 0;

      // Current_node and adjacent node are not connected, so no match.
      *pcount = 0;
#ifdef RULE_CHECK_DEBUG
      DEBUG_PRINTF(10, "component_permute_check failed: "
                   "comp:(r:%d i:%d n:%d) adj:(r:%d i:%d n:%d) ",
                   component_rule->index, component_index, component_node,
                   adjacent_rule->index, adjacent_index, adjacent_node);
      if (CUSTOM_OUTPUT_CHECK(10))
      {
        putset(custom_get_fh(OUTPUT_DEBUG), component_row, &curlen, 80,
               specie->cg->m, FALSE);
      }
      DEBUG_PRINTF(10, "\n");
#endif

      // No error occured.
      return TRUE;
    }
  }

  (*pcount)++;

  // All component rules are satisfied.
  return TRUE;
}

/**
 * Find all nodes that match the given molecule rule and the associated
 * components.  The possibility of component-component bonds will be
 * checked, but the absolutely correct permutations will be found in
 * rule_molecule_permute_check.
 * If update is 1, the list of matched nodes will be save in the
 * component and adjacent rules for later permuting.
 * Only call with update=1 when the entire rule has been confirmed to
 * be correct.
 *
 * The parameter *pcount will contain the number of molecules matched.
 */
static int
rule_molecule_check(specie_t *specie, rule_t *molecule_rule, int *pcount,
                    int update)
{
  int molecule_node;
  int passed;
  rule_t *component_rule;
  int status;
  int n;

  *pcount = 0;
  passed = 0;

  // Create perm_set object to collect all component permutations.
  if (update)
  {
    rule_perm_set_create(molecule_rule);
  }

  // Loop through every node in graph to find matching molecules
  n = specie->cg->n;
  for (molecule_node = 0; molecule_node < n; molecule_node++)
  {
    // Check if molecule rule matches.
#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(10, "Checking molecule rule %d (node %d)\n",
                 molecule_rule->index, molecule_node);
#endif

    if (FALSE == rule_attribute_check(molecule_node, specie, molecule_rule))
    {
      // Rule didn't match.
#ifdef RULE_CHECK_DEBUG
      DEBUG_PRINTF(10, "Molecule attribute check failed\n");
#endif
      continue;
    }

    if (NULL == molecule_rule->adjacent_rule)
    {
      // No component rules, so molecule matches.
      passed = 1;
    }
    else
    {
      // Loop through component_rules.
      for (component_rule = molecule_rule->adjacent_rule;
           component_rule != NULL;
           component_rule = component_rule->and_rule)
      {
        int status;
        int count;

        status = rule_component_check(specie, component_rule, molecule_node,
                                      &count, FALSE);
        if (FALSE == status) return FALSE;

        if (FALSE == rule_count_check(component_rule, count))
        {
          // Component rule didn't find correct number of components.
#ifdef RULE_CHECK_DEBUG
          DEBUG_PRINTF(10, "Component rule (%d) count check failed\n",
                       component_rule->index);
#endif
          passed = 0;
          break;
        }
        else
        {
          // Component rule has correct count, which may be 0.
          passed = 1;
        }
      }
    }

    if (passed != 0)
    {
      // All component rules were satisfied.

      *pcount += 1;

      if (1 == update)
      {
        int perm_base;

        // Don't permute off of permutation sets from different molecule nodes.
        perm_base = perm_set_count(molecule_rule->perm_set);

        // Save node maps in all child rules.
        for (component_rule = molecule_rule->adjacent_rule;
             component_rule != NULL;
             component_rule = component_rule->and_rule)
        {
          int status;
          int count;

          status = rule_component_check(specie, component_rule, molecule_node,
                                        &count, TRUE);
          if (FALSE == status) return FALSE;

#ifndef NDEBUG
          if (FALSE == rule_count_check(component_rule, count))
          {
            // This should never happen.
            assert(0);
            return FALSE;
          }
#endif
        }

        // Generate all permutations for component rule to node mappings
        // for current molecule rule.
        status = rule_perm_set_add_molecule(perm_base, molecule_rule,
                                            molecule_node);
        if (FALSE == status)
        {
          // Not enough nodes to permute rules properly, so rule failed.
          *pcount = 0;
#ifdef RULE_CHECK_DEBUG
          DEBUG_PRINTF(10, "Not enough nodes to permute rules properly\n");
#endif

          return TRUE;
        }

        // Reset component node maps.
        for (component_rule = molecule_rule->adjacent_rule;
             component_rule != NULL;
             component_rule = component_rule->and_rule)
        {
          component_rule->node_map->count = 0;
        }
      }
    }
  }

  return TRUE;
}

/**
 * Create perm_set object of the correct size to hold all permutations
 * of assignments of component nodes to component rules.
 */
static int
rule_perm_set_create(rule_t *molecule_rule)
{
  int i;
  rule_t *component_rule;
  perm_set_t *perm_set;
  int set_size;

  // Count number of component rules molecule rule has
  set_size = 1; // Add one intially for molecule entry.
  component_rule = molecule_rule->adjacent_rule;
  while(component_rule != NULL)
  {
    set_size++;
    component_rule = component_rule->and_rule;
  }

  perm_set = perm_set_create(set_size);

  molecule_rule->perm_set = perm_set;
  molecule_rule->perm_set_index = 0;

  // Loop through component rules.
  // Index 0 will be used by molecule.
  for (i = 1, component_rule = molecule_rule->adjacent_rule;
       component_rule != NULL;
       i++, component_rule = component_rule->and_rule)
  {
    component_rule->perm_set = perm_set;
    component_rule->perm_set_index = i;
  }

  return 1;
}

/**
 * Add component permutations for current molecule.
 */
static int
rule_perm_set_add_molecule(const int perm_base,
                           rule_t *molecule_rule, const int molecule_node)
{
  rule_t *component_rule;
  int status;

  if (NULL == molecule_rule->adjacent_rule)
  {
    // No components present, so simply add entry for molecule node.
    status = perm_set_add_dynints(perm_base, molecule_node,
                                  0, NULL, molecule_rule->perm_set);
    if (0 == status) return 0;

    return 1;
  }

  // Loop through component rules.
  // Index 0 is already in use by molecule.
  for (component_rule = molecule_rule->adjacent_rule;
       component_rule != NULL;
       component_rule = component_rule->and_rule)
  {
    status = perm_set_add_dynints(perm_base, molecule_node,
                                  component_rule->perm_set_index,
                                  component_rule->node_map,
                                  component_rule->perm_set);
    if (0 == status) return 0;
  }

  return 1;
}

/**
 * Find all nodes that match the given component rule.
 * The number of nodes matched is returned in *pcount.
 * If update is 1, the list of matched component nodes are
 * added to rule->node_map.
 * Only call with update as 1 when the parent molecules and other
 * components are known to be good.
 */
static int
rule_component_check(specie_t *specie, rule_t *component_rule,
                     const int molecule_node, int *pcount, const int update)
{
  int component_node;
  set *component_row;
  int component_connection_node;
  set *component_connection_row;
  int passed;

  *pcount = 0;

  component_row = GRAPHROW(specie->cg->g.a, molecule_node, specie->cg->m);

  // Loop through each adjacent node and attempt to match component rule.
  for (component_node = -1;
       (component_node =
          nextelement(component_row, specie->cg->m, component_node)) >= 0;
      )
  {
#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(10, "Checking component rule %d (m_node:%d, c_node:%d)\n",
                 component_rule->index, molecule_node, component_node);
#endif

    if (FALSE == rule_attribute_check(component_node, specie, component_rule))
    {
      // Initial check failed.
#ifdef RULE_CHECK_DEBUG
      DEBUG_PRINTF(10, "Component attribute check failed\n");
#endif
      continue;
    }

    component_connection_row =
      GRAPHROW(specie->cg->g.a, component_node, specie->cg->m);

    // Don't check component-component connections if no rules exist.
    if (NULL == component_rule->adjacent_rule)
    {
      passed = 0;

      // No connections specified, so no connections should be present.
      for (component_connection_node = -1;
           (component_connection_node =
              nextelement(component_connection_row, specie->cg->m,
                          component_connection_node)) >= 0;
          )
      {
        if (component_connection_node != molecule_node)
        {
          // Found component-component connection, so current node fails rule.
          // Proceed to next component.
          passed = 1;
          break;
        } // If component_connection_node != molecule_node.
      } // For each component_connection_node.

      if (1 == passed)
      {
        // Found component-component connection when there shouldn't be any.
#ifdef RULE_CHECK_DEBUG
        DEBUG_PRINTF(10, "Found component-component connection (%d:%d) so "
                     "component fails check\n", component_node,
                     component_connection_node);
#endif
        passed = 0;
      }
      else
      {
        // No connections found, so component was matched.
        passed = 1;
      }
    } // If NULL == component_rule->adjacent_rule.
    else
    {
      int status;
      int count;

      // Component_connection_rule->and_rule is not valid for this molecule.
      // Only one component connection can be preset for each component.
      status = rule_component_connection_check
                 (specie, molecule_node, component_rule->adjacent_rule,
                  component_connection_row, &count, FALSE);
      if (FALSE == status) return FALSE;

      if (FALSE == rule_count_check(component_rule->adjacent_rule, count))
      {
        // Component-component connection rule check failed for current node.
        // Proceed to next node.
#ifdef RULE_CHECK_DEBUG
        DEBUG_PRINTF(10, "Component_connection_check failed\n");
#endif
        continue;
      }
      else
      {
        passed = 1;
      }
    } // If NULL != component_rule->adjacent_node.
                             
    if (1 == passed)
    {
      // This component rule has been satisfied.
      (*pcount)++;

      if (1 == update)
      {
        int status;

        // Parent rules are completely satisfied, so safe to update rules.

        // Add node to component_rule.
        status = dynints_int_add(component_node, component_rule->node_map);
        if (FALSE == status) return FALSE;
      } // If 1 == update.
    } // If count != 0.
  } // For each component node.

  return TRUE;
}

/**
 * Find all nodes that match the given component connection rule.
 * The number of nodes matched is returned in *pcount.
 * If update=1, the list of matched component nodes are
 * added to rule->node_map, 
 * On the first call to this function, update should be 0.
 * Only call with update=1 when the parent molecules and components
 * are known to be good.
 */
static int
rule_component_connection_check(specie_t *specie, const int molecule_node,
                                rule_t *component_connection_rule,
                                set *component_connection_row, int *pcount,
                                const int update)
{
  int component_connection_node;

  *pcount = 0;

  // Loop through each adjacent node and attempt to
  // match component connection rule.
  for (component_connection_node = -1;
       (component_connection_node =
          nextelement(component_connection_row, specie->cg->m,
                      component_connection_node)) >= 0;
      )
  {
#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(10, "Checking component connection rule "
                 "%d (m_node:%d, cc_node:%d)\n",
                 component_connection_rule->index,
                 molecule_node, component_connection_node);
#endif
    // Check if rule matches.
    if ((component_connection_node != molecule_node) &&
        (TRUE == rule_attribute_check(component_connection_node,
                                     specie, component_connection_rule)))
    {
      // Found something, so this rule has been initially satisfied.
      (*pcount)++;
    }
  }

  return TRUE;
}

/**
 * A component is considered matching the rule when the color and all
 * states are present.
 */
static int
rule_attribute_check(int node, specie_t *specie, rule_t *rule)
{
  cgraph_t *cg;

  cg = specie->cg;

  if (((rule->color_value == COLOR_WILDCARD) ||
       (cg->colors.a[node] == rule->color_value)) &&
      ((rule->state_dynints->count == 0) ||
       (1 == dynints_compare_cond(rule->state_dynints,
                                  cg->state_dynints_array.a[node]))))
  {
    return TRUE;
  }

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(12,
               "Attribute check failed: node:%d, specie:%d, rule:%d, "
               "rule_color:%d, node_color:%d "
               "rule_state_count:%ld node_state_count:%ld\n",
               node, specie->id + 1,
               rule->index, rule->color_value, cg->colors.a[node],
               rule->state_dynints->count,
               cg->state_dynints_array.a[node]->count);
#endif
  return FALSE;
}

boolean
rule_setup(const int index, const int color_value, int modifier_flags,
           const count_flag_t count_flag, const int count,
           const double delta, const int xor_index,
           const int or_index, const int adjacent_index, const int and_index,
           rule_t *new_rule)
{
  // Majority of members are set to NULL, so initialize accordingly.
  memset(new_rule, 0, sizeof(rule_t));

  // Set structure members.
  new_rule->index = index;
  new_rule->color_value = color_value;
  new_rule->modifier_flags = modifier_flags;
  new_rule->count_flag = count_flag;
  new_rule->count = count;
  //new_rule->delta = delta;

  new_rule->xor_index = xor_index;
  new_rule->or_index = or_index;
  new_rule->adjacent_index = adjacent_index;
  new_rule->and_index = and_index;
  new_rule->node_map = dynints_create(INITIAL_DYNINTS);
  if (NULL == new_rule->node_map)
  {
    return FALSE;
  }

  new_rule->state_dynints = dynints_create(INITIAL_DYNINTS);
  if (NULL == new_rule->state_dynints)
  {
    return FALSE;
  }

  return TRUE;
}

int
rule_add_bind_index(const int index, rule_t *rule)
{
  rule->attributes_modified = TRUE;

  return ilist_add(index, &rule->bind_index_list);
}

int
rule_add_primary_self_bind_index(const int index, rule_t *rule)
{
  rule->attributes_modified = TRUE;

  return ilist_add(index, &rule->primary_self_bind_index_list);
}

int
rule_add_secondary_self_bind_index(const int index, rule_t *rule)
{
  rule->attributes_modified = TRUE;

  return ilist_add(index, &rule->secondary_self_bind_index_list);
}

int
rule_add_primary_split_index(const int index, rule_t *rule)
{
  rule->attributes_modified = TRUE;

  return ilist_add(index, &rule->primary_split_index_list);
}

int
rule_add_secondary_split_index(const int index, rule_t *rule)
{
  rule->attributes_modified = TRUE;

  return ilist_add(index, &rule->secondary_split_index_list);
}

int
rule_add_state(const int state, const int new_state, rule_t *rule)
{
  int status;

  status = dynints_int_add(state, rule->state_dynints);
  if (FALSE == status)
  {
    return FALSE;
  }

  return status;
}

int
rule_add_modify_state(const int index, const int state,
                      const int new_state, rule_t *rule)
{
  int status;

  rule->attributes_modified = TRUE;

  status = ilist_add(state, &rule->state_change_old_list);
  if (FALSE == status)
  {
    return FALSE;
  }

  status = ilist_add(new_state, &rule->state_change_new_list);
  if (FALSE == status)
  {
    return FALSE;
  }

  return ilist_add(index, &rule->state_change_index_list);
}

void
rule_output(const output_type ptype, const int level, rule_t *rule,
            dyncolors_t *dyncolors, const char *prefix, boolean *used)
{
  int status;
  rule_t *current_rule;

  status = rule_output_raw(ptype, level, rule, dyncolors, prefix, used);

  for (current_rule = rule->or_rule;
       NULL != current_rule;
       current_rule = current_rule->or_rule)
  {
    CUSTOM_PRINTF(ptype, level, "%sOR:\n", prefix);
    status = rule_output_raw(ptype, level, current_rule, dyncolors,
                             prefix, used);
    if (FALSE == status)
    {
      // Ran into rule that was already output, so do not continue.
      break;
    }
  }

}

static int
rule_output_raw(const output_type ptype, const int level, rule_t *rule,
                dyncolors_t *dyncolors, const char *prefix, boolean *used)
{
  char *count_flag_string;
  size_t prefix_len;
  char *new_prefix;
  rule_t *current_rule;
  size_t count;
  int i;
  int status;

  switch(rule->count_flag)
  {
    case LESS_EQ:
      count_flag_string = "<=";
      break;
    case EQUAL:
      count_flag_string = "=";
      break;
    case GREAT_EQ:
      count_flag_string = ">=";
      break;
    default:
      count_flag_string = "";
      break;
  }

  prefix_len = strlen(prefix);

  new_prefix = malloc(prefix_len + 6 + 1);

  memcpy(new_prefix, prefix, prefix_len);
  new_prefix[prefix_len + 0] = ' ';
  new_prefix[prefix_len + 1] = ' ';
  new_prefix[prefix_len + 2] = ' ';
  new_prefix[prefix_len + 3] = ' ';
  new_prefix[prefix_len + 4] = ' ';
  new_prefix[prefix_len + 5] = ' ';
  new_prefix[prefix_len + 6] = '\0';

  if (NULL == rule->xor_rule)
  {
    // No XOR rule, so don't need associated prefix.
    new_prefix[strlen(new_prefix) - 2] = '\0';
  }

  if (NULL == rule->and_rule)
  {
    // No AND rule, so don't need associated prefix.
    new_prefix[strlen(new_prefix) - 2] = '\0';
  }

  new_prefix[strlen(new_prefix) - 2] = '\0';
  if (TRUE == used[rule->index])
  {
    CUSTOM_PRINTF(ptype, level, "%sRule %d: USED\n", new_prefix, rule->index);
    free(new_prefix);
    return FALSE;
  }
  else
  {
    used[rule->index] = TRUE;
  }
  CUSTOM_PRINTF
         (ptype, level, "%sRule %d: 0x%.2X (%s)%s%d (",
          new_prefix, rule->index,
          rule->modifier_flags,
          dyncolors_get_name_from_value(rule->color_value, dyncolors),
          count_flag_string, rule->count);

  count = rule->state_dynints->count;
  if (count > 0)
  {
    CUSTOM_PRINTF(ptype, level, "st:[ ");
    for (i = 0; i < count; i++)
    {
      int value;
      ilist_t *index_list;
      ilist_t *old_list;
      ilist_t *new_list;
      int changed;

      value = rule->state_dynints->ints.a[i];

      changed = 0;
      for (index_list = rule->state_change_index_list,
           old_list = rule->state_change_old_list,
           new_list = rule->state_change_new_list;
           (index_list != NULL) && (old_list != NULL) && (new_list != NULL);
           index_list = index_list->next,
           old_list = old_list->next,
           new_list = new_list->next)
      {
        if (value == old_list->index)
        {
          changed = 1;
          break;
        }
      }


      if (0 == changed)
      {
        // State is not changed.
        CUSTOM_PRINTF(ptype, level, "%s ",
                      dyncolors_get_name_from_value(value, dyncolors));
      }
      else
      {
        // State is changed.
        CUSTOM_PRINTF
          (ptype, level, "%d:%s->%s ", index_list->index,
           dyncolors_get_name_from_value(old_list->index, dyncolors),
           dyncolors_get_name_from_value(new_list->index, dyncolors));
      }
    }
    CUSTOM_PRINTF(ptype, level, "]");
  }

  rule_ilist_output(ptype, level,
                         "ps", rule->primary_split_index_list);
  rule_ilist_output(ptype, level,
                         "ss", rule->secondary_split_index_list);
  rule_ilist_output(ptype, level,
                         "b", rule->bind_index_list);
  rule_ilist_output(ptype, level,
                         "psb", rule->primary_self_bind_index_list);
  rule_ilist_output(ptype, level,
                         "ssb", rule->secondary_self_bind_index_list);
  rule_ilist_output(ptype, level,
                         "tb", rule->transient_bind_index_list);
  rule_ilist_output(ptype, level, "tpsb",
                         rule->primary_transient_self_bind_index_list);
  rule_ilist_output(ptype, level, "tssb",
                         rule->secondary_transient_self_bind_index_list);

  CUSTOM_PRINTF(ptype, level, ") { ");
  if (-1 != rule->xor_index)
  {
    CUSTOM_PRINTF(ptype, level, "xor:%d ", rule->xor_index);
  }
  if (-1 != rule->or_index)
  {
    CUSTOM_PRINTF(ptype, level, "or:%d ", rule->or_index);
  }
  if (-1 != rule->adjacent_index)
  {
    CUSTOM_PRINTF(ptype, level, "adj:%d ", rule->adjacent_index);
  }
  if (-1 != rule->and_index)
  {
    CUSTOM_PRINTF(ptype, level, "and:%d ", rule->and_index);
  }
  CUSTOM_PRINTF(ptype, level, "} <");

  {
    color_type_t color_type;

    color_type = dyncolors_get_type_from_value(rule->color_value, dyncolors);
    if (COLOR_TYPE_MOLECULE == color_type)
    {
      // Output all component permutations.
      if (NULL == rule->perm_set)
      {
        CUSTOM_PRINTF(ptype, level, "perm_set:NULL");
      }
      else
      {
        if (NULL == rule->perm_set->current_perm)
        {
          CUSTOM_PRINTF(ptype, level, "index:%d node:NULL perm_set:",
                        rule->perm_set_index);
        }
        else
        {
          CUSTOM_PRINTF(ptype, level, "index:%d node:%d perm_set:",
                        rule->perm_set_index,
                        rule->perm_set->current_perm[rule->perm_set_index]);
        }

        PERM_SET_OUTPUT(ptype, level, rule->perm_set);
      }
    }
    else if (COLOR_TYPE_COMPONENT == color_type)
    {
      int i;
      int *pnode;
      int map_count;

      // Output index of component in permutations.
      CUSTOM_PRINTF(ptype, level, "index:%d map:[", rule->perm_set_index);

      // Output node map.
      map_count = rule->node_map->count;
      for (i = 0, pnode = rule->node_map->ints.a;
           i < map_count;
           i++, pnode++)
      {
        if (i != 0)
        {
          CUSTOM_PRINTF(ptype, level, ",");
        }
        CUSTOM_PRINTF(ptype, level, "%d", *pnode);
      }
      CUSTOM_PRINTF(ptype, level, "]");
    }
  }

  CUSTOM_PRINTF(ptype, level, ">\n");

  new_prefix[strlen(new_prefix)] = ' ';

  if (NULL != rule->adjacent_rule)
  {
    CUSTOM_PRINTF(ptype, level, "%sAdjacent:\n", new_prefix);
    for (current_rule = rule->adjacent_rule;
         NULL != current_rule;
         current_rule = current_rule->and_rule)
    {
      status = rule_output_raw(ptype, level,current_rule,
                               dyncolors, new_prefix, used);
      if (FALSE == status)
      {
        // Ran into rule that was already output, so do not continue.
        break;
      }
    }
  }

  new_prefix[strlen(new_prefix) - 2] = '\0';

  if (NULL != rule->and_rule)
  {
    for (current_rule = rule->and_rule;
         NULL != current_rule;
         current_rule = current_rule->and_rule)
    {
      new_prefix[strlen(new_prefix) - 2] = '\0';
      CUSTOM_PRINTF(ptype, level, "%sAND:\n", new_prefix);
      new_prefix[strlen(new_prefix)] = ' ';

      status = rule_output_raw(ptype, level, current_rule,
                               dyncolors, new_prefix, used);
      if (FALSE == status)
      {
        // Ran into rule that was already output, so do not continue.
        break;
      }
    }

    new_prefix[strlen(new_prefix) - 2] = '\0';
  }

  if (NULL != rule->xor_rule)
  {
    CUSTOM_PRINTF(ptype, level, "%sXOR:\n", prefix);
    for (current_rule = rule->xor_rule;
         NULL != current_rule;
         current_rule = current_rule->and_rule)
    {
      new_prefix[strlen(new_prefix) - 2] = '\0';
      CUSTOM_PRINTF(ptype, level, "%sXOR:\n", new_prefix);
      new_prefix[strlen(new_prefix)] = ' ';

      status = rule_output_raw(ptype, level, current_rule, dyncolors,
                               new_prefix, used);
      if (FALSE == status)
      {
        // Ran into rule that was already output, so do not continue.
        break;
      }
    }

    new_prefix[strlen(new_prefix) - 2] = '\0';
  }

  free (new_prefix);

  return TRUE;
}

int
rule_get_max_primary_split_elements(rule_t *rule, char *used_rules)
{
  int elements;
  int other_elements;
  ilist_t *current_ilist;

  if ((NULL == rule) || (TRUE == used_rules[rule->index]))
  {
    return 0;
  }

  // Mark rule already checked.
  used_rules[rule->index] = TRUE;

  elements = 0;

  current_ilist = rule->primary_split_index_list;
  while (NULL != current_ilist)
  {
    if ((current_ilist->index + 1) > elements)
    {
      elements = current_ilist->index + 1;
    }

    current_ilist = current_ilist->next;
  }

  other_elements = rule_get_max_primary_split_elements(rule->xor_rule,
                                                       used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_primary_split_elements(rule->or_rule,
                                                       used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_primary_split_elements(rule->adjacent_rule,
                                                       used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_primary_split_elements(rule->and_rule,
                                                       used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  return elements;
}

int
rule_get_max_secondary_split_elements(rule_t *rule, char *used_rules)
{
  int elements;
  int other_elements;
  ilist_t *current_ilist;

  if ((NULL == rule) || (TRUE == used_rules[rule->index]))
  {
    return 0;
  }

  // Mark rule already checked.
  used_rules[rule->index] = TRUE;

  elements = 0;

  current_ilist = rule->secondary_split_index_list;
  while (NULL != current_ilist)
  {
    if ((current_ilist->index + 1) > elements)
    {
      elements = current_ilist->index + 1;
    }

    current_ilist = current_ilist->next;
  }

  other_elements = rule_get_max_secondary_split_elements(rule->xor_rule,
                                                         used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_secondary_split_elements(rule->or_rule,
                                                         used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_secondary_split_elements(rule->adjacent_rule,
                                                         used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_secondary_split_elements(rule->and_rule,
                                                         used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  return elements;
}

int
rule_get_max_bind_elements(rule_t *rule, char *used_rules)
{
  int elements;
  int other_elements;
  ilist_t *current_ilist;

  if ((NULL == rule) || (TRUE == used_rules[rule->index]))
  {
    return 0;
  }

  // Mark rule already checked.
  used_rules[rule->index] = TRUE;

  elements = 0;

  current_ilist = rule->bind_index_list;
  while (NULL != current_ilist)
  {
    if ((current_ilist->index + 1) > elements)
    {
      elements = current_ilist->index + 1;
    }

    current_ilist = current_ilist->next;
  }

  other_elements = rule_get_max_bind_elements(rule->xor_rule,
                                              used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_bind_elements(rule->or_rule,
                                              used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_bind_elements(rule->adjacent_rule,
                                              used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_bind_elements(rule->and_rule,
                                              used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  return elements;
}

int
rule_get_max_primary_self_bind_elements(rule_t *rule, char *used_rules)
{
  int elements;
  int other_elements;
  ilist_t *current_ilist;

  if ((NULL == rule) || (TRUE == used_rules[rule->index]))
  {
    return 0;
  }

  // Mark rule already checked.
  used_rules[rule->index] = TRUE;

  elements = 0;

  current_ilist = rule->primary_self_bind_index_list;
  while (NULL != current_ilist)
  {
    if ((current_ilist->index + 1) > elements)
    {
      elements = current_ilist->index + 1;
    }

    current_ilist = current_ilist->next;
  }

  other_elements = rule_get_max_primary_self_bind_elements(rule->xor_rule,
                                                           used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_primary_self_bind_elements(rule->or_rule,
                                                           used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_primary_self_bind_elements(rule->adjacent_rule,
                                                           used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_primary_self_bind_elements(rule->and_rule,
                                                           used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  return elements;
}

int
rule_get_max_secondary_self_bind_elements(rule_t *rule, char *used_rules)
{
  int elements;
  int other_elements;
  ilist_t *current_ilist;

  if ((NULL == rule) || (TRUE == used_rules[rule->index]))
  {
    return 0;
  }

  // Mark rule already checked.
  used_rules[rule->index] = TRUE;

  elements = 0;

  current_ilist = rule->secondary_self_bind_index_list;
  while (NULL != current_ilist)
  {
    if ((current_ilist->index + 1) > elements)
    {
      elements = current_ilist->index + 1;
    }

    current_ilist = current_ilist->next;
  }

  other_elements = rule_get_max_secondary_self_bind_elements(rule->xor_rule,
                                                             used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_secondary_self_bind_elements(rule->or_rule,
                                                             used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_secondary_self_bind_elements
                     (rule->adjacent_rule, used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_secondary_self_bind_elements(rule->and_rule,
                                                             used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  return elements;
}

int
rule_get_max_state_change_elements(rule_t *rule, char *used_rules)
{
  int elements;
  int other_elements;
  ilist_t *current_ilist;

  if ((NULL == rule) || (TRUE == used_rules[rule->index]))
  {
    return 0;
  }

  // Mark rule already checked.
  used_rules[rule->index] = TRUE;

  elements = 0;

  current_ilist = rule->state_change_index_list;
  while (NULL != current_ilist)
  {
    if ((current_ilist->index + 1) > elements)
    {
      elements = current_ilist->index + 1;
    }

    current_ilist = current_ilist->next;
  }

  other_elements = rule_get_max_state_change_elements(rule->xor_rule,
                                                      used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_state_change_elements(rule->or_rule,
                                                      used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_state_change_elements
                     (rule->adjacent_rule, used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_state_change_elements(rule->and_rule,
                                                      used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  return elements;
}

int
rule_get_max_transient_bind_elements(rule_t *rule, char *used_rules)
{
  int elements;
  int other_elements;
  ilist_t *current_ilist;

  if ((NULL == rule) || (TRUE == used_rules[rule->index]))
  {
    return 0;
  }

  // Mark rule already checked.
  used_rules[rule->index] = TRUE;

  elements = 0;

  current_ilist = rule->transient_bind_index_list;
  while (NULL != current_ilist)
  {
    if ((current_ilist->index + 1) > elements)
    {
      elements = current_ilist->index + 1;
    }

    current_ilist = current_ilist->next;
  }

  other_elements = rule_get_max_transient_bind_elements(rule->xor_rule,
                                                        used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_transient_bind_elements(rule->or_rule,
                                                        used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_transient_bind_elements(rule->adjacent_rule,
                                                        used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements = rule_get_max_transient_bind_elements(rule->and_rule,
                                                        used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  return elements;
}

int
rule_get_max_primary_transient_self_bind_elements(rule_t *rule,
                                                  char *used_rules)
{
  int elements;
  int other_elements;
  ilist_t *current_ilist;

  if ((NULL == rule) || (TRUE == used_rules[rule->index]))
  {
    return 0;
  }

  // Mark rule already checked.
  used_rules[rule->index] = TRUE;

  elements = 0;

  current_ilist = rule->primary_transient_self_bind_index_list;
  while (NULL != current_ilist)
  {
    if ((current_ilist->index + 1) > elements)
    {
      elements = current_ilist->index + 1;
    }

    current_ilist = current_ilist->next;
  }

  other_elements =
    rule_get_max_primary_transient_self_bind_elements(rule->xor_rule,
                                                      used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements =
    rule_get_max_primary_transient_self_bind_elements(rule->or_rule,
                                                      used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements =
    rule_get_max_primary_transient_self_bind_elements(rule->adjacent_rule,
                                                      used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements =
    rule_get_max_primary_transient_self_bind_elements(rule->and_rule,
                                                      used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  return elements;
}

int
rule_get_max_secondary_transient_self_bind_elements(rule_t *rule,
                                                    char *used_rules)
{
  int elements;
  int other_elements;
  ilist_t *current_ilist;

  if ((NULL == rule) || (TRUE == used_rules[rule->index]))
  {
    return 0;
  }

  // Mark rule already checked.
  used_rules[rule->index] = TRUE;

  elements = 0;

  current_ilist =
    rule->secondary_transient_self_bind_index_list;
  while (NULL != current_ilist)
  {
    if ((current_ilist->index + 1) > elements)
    {
      elements = current_ilist->index + 1;
    }

    current_ilist = current_ilist->next;
  }

  other_elements =
    rule_get_max_secondary_transient_self_bind_elements(rule->xor_rule,
                                                        used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements =
    rule_get_max_secondary_transient_self_bind_elements(rule->or_rule,
                                                        used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements =
    rule_get_max_secondary_transient_self_bind_elements(rule->adjacent_rule,
                                                        used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  other_elements =
    rule_get_max_secondary_transient_self_bind_elements(rule->and_rule,
                                                        used_rules);
  if (other_elements > elements)
  {
    elements = other_elements;
  }

  return elements;
}

int
rule_create_specie(rule_t *rule, int *pnew_specie_id, dynspecies_t *dynspecies,
                   dynrules_t *dynrules)
{
  int node_count;
  cgraph_t *new_cgraph;
  dynmaps_t *rule_node_maps;
  dynints_t **state_dynints_array;
  rule_t *molecule_rule;
  rule_t *component_rule;
  int molecule_node;
  int component_node;
  int status;
  int new;
  rule_t *rules;

  rules = dynrules->rules.a;

  rule_node_maps = dynmaps_create(INITIAL_DYNMAPS);

  // Build new specie as cgraph.
  new_cgraph = cgraph_create(INITIAL_CGRAPH_SIZE);

  // Count all nodes, map rules to nodes.
  node_count = 0;
  for (molecule_rule = rule;
       NULL != molecule_rule;
       molecule_rule = rules + molecule_rule->and_index)
  {
    // Add node to graph.
    cgraph_add_node(new_cgraph, molecule_rule->color_value);

    // Molecule node.
    node_count++;

    if (molecule_rule->adjacent_index >= 0)
    {
      for (component_rule = rules + molecule_rule->adjacent_index;
           component_rule != NULL;
           component_rule = rules + component_rule->and_index)
      {
        // Save rule - node association.
        dynmaps_map_add(component_rule->index, node_count, rule_node_maps);

        // Add node to graph.
        cgraph_add_node(new_cgraph, component_rule->color_value);

        // Component node.
        node_count++;

        if (component_rule->and_index < 0)
        {
          break;
        }
      }
    }

    if (molecule_rule->and_index < 0)
    {
      break;
    }
  }

  state_dynints_array = new_cgraph->state_dynints_array.a;

  // Create connections and states in cgraph.
  node_count = 0;
  for (molecule_rule = rule;
       NULL != molecule_rule;
       molecule_rule = rules + molecule_rule->and_index)
  {
    molecule_node = node_count;
                 
    // Molecule node.
    node_count++;

    if (molecule_rule->adjacent_index >= 0)
    {
      for (component_rule = rules + molecule_rule->adjacent_index;
           component_rule != NULL;
           component_rule = rules + component_rule->and_index)
      {
        dynints_t *component_states;
        int adj_component;

        component_node = node_count;

        component_states = component_rule->state_dynints;

        if (component_states != NULL)
        {
          // Save state(s).
          status = dynints_copy(component_rule->state_dynints,
                                *(state_dynints_array + component_node));
          if (0 == status)
          {
            dynmaps_destroy(rule_node_maps);
            cgraph_destroy(new_cgraph);
            return FALSE;
          }
        }

        // Component is always adjacent to molecule.
        cgraph_connect_nodes(new_cgraph, molecule_node, component_node);

        // Add any connection requirements.
        adj_component = component_rule->adjacent_index;
        if ((adj_component >= 0) && (adj_component < component_rule->index))
        {
          // Add connection between components.
          cgraph_connect_nodes(new_cgraph, adj_component,
                               component_rule->index);
        }

        // Component node.
        node_count++;

        if (component_rule->and_index < 0)
        {
          break;
        }
      }
    }

    if (molecule_rule->and_index < 0)
    {
      break;
    }
  }

  /* Normalize color partitions. */
  status = cgraph_normalize_partitions(new_cgraph);
  if (0 == status)
  {
    return 0;
  }

  /* Add specie */
  status = dynspecies_specie_add(new_cgraph, "created",
                                 0, TRUE, TRUE, FALSE, TRUE,
                                 pnew_specie_id, &new, world->dynspecies,
                                 world);
  if (0 == status)
  {
    dynmaps_destroy(rule_node_maps);
    cgraph_destroy(new_cgraph);
    return FALSE;
  }

  dynmaps_destroy(rule_node_maps);
  cgraph_destroy(new_cgraph);

  return 1;
}

int
rule_complete(rule_t *rule, dynrules_t *dynrules)
{
  if (-1 != rule->xor_index)
  {
    rule->xor_rule = dynrules->rules.a + rule->xor_index;
  }
  else
  {
    rule->xor_rule = NULL;
  }

  if (-1 != rule->or_index)
  {
    rule->or_rule = dynrules->rules.a + rule->or_index;
  }
  else
  {
    rule->or_rule = NULL;
  }

  if (-1 != rule->adjacent_index)
  {
    rule->adjacent_rule = dynrules->rules.a + rule->adjacent_index;
  }
  else
  {
    rule->adjacent_rule = NULL;
  }

  if (-1 != rule->and_index)
  {
    rule->and_rule = dynrules->rules.a + rule->and_index;
  }
  else
  {
    rule->and_rule = NULL;
  }

  return 1;
}

static void
rule_ilist_output(const output_type ptype, const int level,
                       const char *label, ilist_t *ilist)
{
  if (NULL != ilist)
  {
    ilist_t *current_ni;

    current_ni = ilist;

    CUSTOM_PRINTF(ptype, level, "%s:[", label);
    for (;;)
    {
      CUSTOM_PRINTF(ptype, level, "%d", current_ni->index);

      current_ni = current_ni->next;
      if (NULL != current_ni)
      {
        CUSTOM_PRINTF(ptype, level, ",");
      }
      else
      {
        break;
      }
    }
    CUSTOM_PRINTF(ptype, level, "]");
  }
}

void
rule_destroy(rule_t *rule)
{
  ilist_destroy(rule->primary_split_index_list);
  ilist_destroy(rule->secondary_split_index_list);
  ilist_destroy(rule->bind_index_list);
  ilist_destroy(rule->primary_self_bind_index_list);
  ilist_destroy(rule->secondary_self_bind_index_list);
  ilist_destroy(rule->state_change_index_list);
  ilist_destroy(rule->state_change_old_list);
  ilist_destroy(rule->state_change_new_list);
  ilist_destroy(rule->transient_bind_index_list);
  ilist_destroy(rule->primary_transient_self_bind_index_list);
  ilist_destroy(rule->secondary_transient_self_bind_index_list);
  dynints_destroy(rule->node_map);

  // Perm_set gets destroyed after every rule check.

  free(rule);
}

/*
** Verify that the current number of elements matched satisfies the rule.
*/
static boolean
rule_count_check(rule_t *rule, const int count)
{
  if (((EQUAL  == rule->count_flag) &&
       (count != rule->count)) ||
      ((LESS_EQ == rule->count_flag) &&
       (count >  rule->count)) ||
      ((GREAT_EQ == rule->count_flag) &&
       (count  <  rule->count)))
  {
#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(10,
               "Failed: Count:%d Rule Count:%d Flag:%.2x.\n", count,
               rule->count, rule->count_flag);
#endif

    return FALSE;
  }
  else
  {
#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(10,
               "Success: rule:%d, Count:%d Rule Count:%d, Flag:%.2x.\n",
               rule->index, count, rule->count, rule->count_flag);
#endif

    return TRUE;
  }
}

