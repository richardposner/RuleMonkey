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
#include "dynarray.h"
#include "constants.h"
#include "output.h"

void suppress_warnings13() {SUPPRESS_WARNINGS;}

// Local function declarations.
static boolean dynrules_resize(dynrules_t *dynrules, const size_t count);
static int dynrules_get_bind_colors(rule_t *rule, dyncolors_t *dyncolors,
                                    int *pmolecule_color_value,
                                    int *pcomponent_color_value);

int
dynrules_same_components(rule_t *rule_a, rule_t *rule_b,
                         dyncolors_t *dyncolors)
{
  int a_molecule;
  int a_component;
  int b_molecule;
  int b_component;
  int status;

  // Get component and molecule colors for first reactant.
  status = dynrules_get_bind_colors(rule_a, dyncolors,
                                    &a_molecule, &a_component);
  if (0 == status)
  {
    return -1;
  }

  // Get component and molecule colors for second reactant.
  status = dynrules_get_bind_colors(rule_b, dyncolors,
                                    &b_molecule, &b_component);
  if (0 == status)
  {
    return -1;
  }

  if ((a_molecule == b_molecule) && (a_component == b_component))
  {
    return 1;
  }

  return 0;
}

static int
dynrules_get_bind_colors(rule_t *rule, dyncolors_t *dyncolors,
                         int *pmolecule_color_value,
                         int *pcomponent_color_value)
{
  rule_t *molecule_rule;
  rule_t *component_rule;

  if (NULL == rule)
  {
    return 0;
  }

  // Cycle through molecules
  for (molecule_rule = rule;
       molecule_rule != NULL;
       molecule_rule = molecule_rule->and_rule)
  {
    // Cycle through components
    for (component_rule = molecule_rule->adjacent_rule;
         NULL != component_rule;
         component_rule = component_rule->and_rule)
    {
      if (component_rule->bind_index_list != NULL)
      {
        // Found bind point.
        *pmolecule_color_value  = molecule_rule->color_value;
        *pcomponent_color_value  = component_rule->color_value;
        return 1;
      }
    }
  }

  // Bind rule not found.
  return 0;
}

boolean
dynrules_rule_create(const int color_value,
                     int modifier_flags,
                     const count_flag_t count_flag, const int count,
                     const double delta,
                     const int xor_index, const int or_index,
                     const int adjacent_index, const int and_index,
                     dynrules_t *dynrules)
{
  rule_t *new_rule;
  int status;

  new_rule = dynrules_rule_alloc(dynrules);
  if (NULL == new_rule)
  {
    return FALSE;
  }

  status = rule_setup((int)dynrules->count - 1, color_value, modifier_flags,
                      count_flag, count, delta, xor_index, or_index,
                      adjacent_index, and_index, new_rule);

  return status;
}

boolean
dynrules_rule_add_modifier_flags(const int index, const int or_mask,
                                 dynrules_t *dynrules)
{
  // Make sure rule exists.
  if (index >= dynrules->count)
  {
    error_printf("Rule with index %d doesn't exist.\n", index);
    return FALSE;
  }

  dynrules->rules.a[index].modifier_flags |= or_mask;

  return TRUE;
}

boolean
dynrules_rule_remove_modifier_flags(const int index,
                                    const int and_mask,
                                    dynrules_t *dynrules)
{
  // Make sure rule exists.
  if (index >= dynrules->count)
  {
    error_printf("Rule with index %d doesn't exist.\n", index);
    return FALSE;
  }

  dynrules->rules.a[index].modifier_flags &= and_mask;

  return TRUE;
}

boolean
dynrules_rule_modify_count(const int index,
                           const count_flag_t count_flag, const int count,
                           dynrules_t *dynrules)
{
  // Make sure rule exists.
  if (index >= dynrules->count)
  {
    error_printf("Rule with index %d doesn't exist.\n", index);
    return FALSE;
  }

  dynrules->rules.a[index].count_flag = count_flag;
  dynrules->rules.a[index].count = count;

  return TRUE;
}

boolean
dynrules_rule_set_reset(const int index, const int value,
                        dynints_t *rules_used, dynrules_t *dynrules)
{
  boolean status;

  if (rules_used->count <= index)
  {
    status = dynints_resize(rules_used, dynrules->count);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }
    rules_used->count = dynrules->count;
  }
  rules_used->ints.a[index] = value;

  return TRUE;
}

boolean
dynrules_rule_connect_zero_or_one_rule(const int index, dynrules_t *dynrules)
{
  int status;
  int new_index;

  status = dynrules_rule_create(COLOR_WILDCARD, NO_MODIFIER_FLAGS, LESS_EQ,
                                1, 1.0, -1, -1, -1, -1, dynrules);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

  new_index = dynrules->count - 1;

  dynrules_rule_add_adjacent_rule(index, new_index, dynrules);

  return TRUE;
}

boolean
dynrules_rule_connect_any_one_rule(const int index, dynrules_t *dynrules)
{
  int status;
  int new_index;

  status = dynrules_rule_create(COLOR_WILDCARD, NO_MODIFIER_FLAGS, EQUAL,
                                1, 1.0, -1, -1, -1, -1, dynrules);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

  new_index = dynrules->count - 1;

  dynrules_rule_add_adjacent_rule(index, new_index, dynrules);

  return TRUE;
}

void
dynrules_rule_add_xor_rule(const int index, const int xor_index,
                           dynrules_t *dynrules)
{
  rule_t *current_rule;

  current_rule = dynrules->rules.a + index;

  if (-1 == current_rule->xor_index)
  {
#ifdef RULE_CREATE_DEBUG
    DEBUG_PRINTF
      (4,
       "add_xor_rule.  index=%d, xor_index=%d, current=%d\n",
       index, xor_index, current_rule->index);
#endif
    current_rule->xor_index = xor_index;

    // Done.
    return;
  }

  // Current rule already has an xor rule, need to get 'and' to
  // this xor rule.
#ifdef RULE_CREATE_DEBUG
  DEBUG_PRINTF(4, "xor:");
#endif
  dynrules_rule_add_and_rule(current_rule->xor_index,
                             xor_index, dynrules);

  return;
}

void
dynrules_rule_add_or_rule(const int index, const int or_index,
                          dynrules_t *dynrules)
{
  rule_t *current_rule;

  current_rule = dynrules->rules.a + index;

  if (-1 == current_rule->or_index)
  {
#ifdef RULE_CREATE_DEBUG
    DEBUG_PRINTF(4,
                 "add_or_rule.  index=%d, or_index=%d, current=%d\n",
                 index, or_index, current_rule->index);
#endif
    current_rule->or_index = or_index;

    // Done.
    return;
  }

  // Current rule already has an or rule, need to get 'and' to
  // this or rule.
#ifdef RULE_CREATE_DEBUG
  DEBUG_PRINTF(4, "or:");
#endif
  dynrules_rule_add_and_rule(current_rule->or_index,
                             or_index, dynrules);

  return;
}

void
dynrules_rule_add_and_rule(const int index, const int and_index,
                           dynrules_t *dynrules)
{
  rule_t *current_rule;

  current_rule = dynrules->rules.a + index;

  while (-1 != current_rule->and_index)
  {
    // Get to the end of the 'and' chain.
    current_rule = dynrules->rules.a + current_rule->and_index;
  }

#ifdef RULE_CREATE_DEBUG
  DEBUG_PRINTF(4,
               "add_and_rule.  index=%d, and_index=%d, current=%d\n",
               index, and_index, current_rule->index);
#endif

  // Add 'and' rule to end of 'and' chain.
  current_rule->and_index = and_index;

  return;
}

void
dynrules_rule_add_adjacent_rule(const int index, const int adjacent_index,
                                dynrules_t *dynrules)
{
  rule_t *current_rule;

  current_rule = dynrules->rules.a + index;

  if (-1 == current_rule->adjacent_index)
  {
#ifdef RULE_CREATE_DEBUG
    DEBUG_PRINTF
      (4,
       "add_adjacent_rule.  index=%d, adjacent_index=%d, current=%d\n",
       index, adjacent_index, current_rule->index);
#endif
    current_rule->adjacent_index = adjacent_index;

    // Done.
    return;
  }

  // Current rule already has an adjacent rule, need to get 'and' to
  // this adjacent rule.
#ifdef RULE_CREATE_DEBUG
  DEBUG_PRINTF(4, "adjacent:");
#endif
  dynrules_rule_add_and_rule(current_rule->adjacent_index,
                              adjacent_index, dynrules);

  return;
}

boolean
dynrules_rule_remove_and_rule(const int index, const int and_index,
                               dynrules_t *dynrules)
{
  rule_t *rules;
  rule_t *current_rule;

  rules = dynrules->rules.a;

  for (current_rule = rules + index;
       -1 != current_rule->and_index;
       current_rule = rules + current_rule->and_index)
  {
    if (and_index == current_rule->and_index)
    {
#ifdef RULE_CREATE_DEBUG
      DEBUG_PRINTF
        (4,
         "remove_and_rule.  index=%d, and_index=%d, current=%d\n",
         index, and_index, current_rule->index);
#endif
      // Found rule to remove.
      current_rule->and_index = rules[and_index].and_index;

      rules[and_index].and_index = -1;

      return TRUE;
    }
  }

  // Fell off end without finding rule to remove.
  assert(0);
  return FALSE;
}

boolean
dynrules_rule_remove_adjacent_rule(const int index, const int adjacent_index,
                                   dynrules_t *dynrules)
{
  rule_t *current_rule;
  int status;

  current_rule = dynrules->rules.a + index;

  if (-1 == current_rule->adjacent_index)
  {
    // Adjacent rule not found.
    assert(0);
    return FALSE;
  }

  if (adjacent_index == current_rule->adjacent_index)
  {
#ifdef RULE_CREATE_DEBUG
    DEBUG_PRINTF(4, "remove_adjacent_rule: "
                 "adjacent_index=%d, current=%d, adjacent_and=%d\n",
                 adjacent_index, current_rule->index,
                 dynrules->rules.a[adjacent_index].and_index);
#endif
    // Found rule to remove.
    current_rule->adjacent_index = dynrules->rules.a[adjacent_index].and_index;

    dynrules->rules.a[adjacent_index].and_index = -1;

    return TRUE;
  }

  // Search for rule.
#ifdef RULE_CREATE_DEBUG
  DEBUG_PRINTF(4, "adjacent, index:%d:", index);
#endif
  status = dynrules_rule_remove_and_rule(current_rule->adjacent_index,
                                          adjacent_index, dynrules);

  return status;
}

rule_t *
dynrules_rule_alloc(dynrules_t *dynrules)
{
  int status;

  status = dynrules_resize(dynrules, dynrules->count + 1);
  if (0 == status)
  {
    assert(0);
    return NULL;
  }

  return &dynrules->rules.a[dynrules->count++];
}

boolean
dynrules_rule_add_bind_index(const int index, const int bind_index,
                             dynrules_t *dynrules)
{
  rule_add_bind_index(bind_index, dynrules->rules.a + index);

  return TRUE;
}

boolean
dynrules_rule_add_primary_self_bind_index(const int index,
                                          const int bind_index,
                                          dynrules_t *dynrules)
{
  rule_add_primary_self_bind_index(bind_index, dynrules->rules.a + index);

  return TRUE;
}

boolean
dynrules_rule_add_secondary_self_bind_index(const int index,
                                            const int bind_index,
                                            dynrules_t *dynrules)
{
  rule_add_secondary_self_bind_index(bind_index, dynrules->rules.a + index);

  return TRUE;
}

boolean
dynrules_rule_add_primary_split_index(const int index, const int bind_index,
                                      dynrules_t *dynrules)
{
  rule_add_primary_split_index(bind_index, dynrules->rules.a + index);

  return TRUE;
}

boolean
dynrules_rule_add_secondary_split_index(const int index,
                                        const int bind_index,
                                        dynrules_t *dynrules)
{
  rule_add_secondary_split_index(bind_index, dynrules->rules.a + index);

  return TRUE;
}

boolean
dynrules_map_molecules(const int lhs_count, const int rhs_count,
                       const int molecule_star_lhs,
                       const int molecule_star_rhs,
                       int *specie_rule_indexes,
                       dynmaps_t *dynmaps,
                       dynints_t *extra_lhs_molecules,
                       dynints_t *extra_rhs_molecules,
                       dynconnections_t *dynconnections,
                       dyncolors_t *dyncolors, dynrules_t *dynrules)
{
  int status;
  int current_rule_index;
  rule_t *current_rule;
  size_t lhs_index;
  size_t rhs_index;
  size_t last_rhs_index;
  int lhs_done;
  int rhs_done;

  lhs_done = 0;
  rhs_done = 0;

  lhs_index = specie_rule_indexes[0];
  rhs_index = specie_rule_indexes[lhs_count];
  last_rhs_index = specie_rule_indexes[lhs_count + rhs_count - 1];

  assert(rhs_index <= dynrules->count);

  current_rule_index = lhs_index;
  current_rule = dynrules->rules.a + lhs_index;

  // Process all lhs rules.
  while (current_rule_index < rhs_index)
  {
    // Find next molecule rule.
    while (COLOR_TYPE_MOLECULE !=
           dyncolors_get_type_from_value(current_rule->color_value,
                                         dyncolors))
    {
      current_rule_index++;
      current_rule++;

      if (current_rule_index >= rhs_index)
      {
        // All lhs molecules have been mapped.
        lhs_done = 1;
        break;
      }
    }

    if ((current_rule_index >= rhs_index) ||
        (1 == lhs_done))
    {
      // All lhs molecules have been mapped.
      break;
    }

    // Found a molecule on the lhs.  Search for matching molecule on rhs.
    status = dynrules_find_matching_molecule(current_rule_index, rhs_index,
                                             dynmaps, dynconnections,
                                             dyncolors, dynrules);
    if (FALSE == status)
    {
      /* Keep track of the molecules that need to be created. */
      status = dynints_int_add(current_rule_index, extra_lhs_molecules);
      if (0 == status)
      {
        return FALSE;
      }
    }

    // Try next lhs molecule.
    current_rule_index++;
    current_rule++;
  }

  // Check that all rhs molecules have been mapped.
  while (current_rule_index <= last_rhs_index)
  {
    // Find next molecule rule.
    while (COLOR_TYPE_MOLECULE !=
           dyncolors_get_type_from_value(current_rule->color_value,
                                         dyncolors))
    {
      current_rule_index++;
      current_rule++;

      if (current_rule_index > last_rhs_index)
      {
        // All molecules have been mapped.
        rhs_done = 1;
        break;
      }
    }

    if ((current_rule_index > last_rhs_index) ||
        (1 == rhs_done))
    {
      // All molecules have been mapped.
      return TRUE;
    }


    // Check if molecule on rhs is mapped to lhs.
    status = dynmaps_dest_exists(current_rule_index, dynmaps);
    if (FALSE == status)
    {
      /* Keep track of the molecules that need to be created. */
      status = dynints_int_add(current_rule_index, extra_rhs_molecules);
      if (0 == status)
      {
        return FALSE;
      }
    }

    // Try next rhs molecule.
    current_rule_index++;
    current_rule++;
  }

  // All molecules have been mapped.
  return TRUE;
}

boolean
dynrules_find_matching_molecule(size_t current_lhs_rule_index,
                                size_t rhs_index, dynmaps_t *dynmaps,
                                dynconnections_t *dynconnections,
                                dyncolors_t *dyncolors, dynrules_t *dynrules)
{
  int status;
  size_t rule_count;
  size_t current_rhs_rule_index;
  rule_t *current_rhs_rule;

  // Loop through each rhs molecule until a match is found.
  rule_count = dynrules->count;
  current_rhs_rule_index = rhs_index;
  current_rhs_rule = dynrules->rules.a + rhs_index;

  // Process all rhs rules.
  while (current_rhs_rule_index < rule_count)
  {
    // Find next molecule rule.
    while (COLOR_TYPE_MOLECULE !=
           dyncolors_get_type_from_value(current_rhs_rule->color_value,
                                         dyncolors))
    {
      current_rhs_rule_index++;
      current_rhs_rule++;

      if (current_rhs_rule_index >= rule_count)
      {
        // The lhs molecule didn't match any of the rhs molecules.
        error_printf("lhs molecule didn't match any rhs molecule, "
                     "ran out of rules to check.\n");
        return FALSE;
      }
    }

    // Found a molecule on the rhs.  Make sure it hasn't been mapped yet.
    status = dynmaps_dest_exists(current_rhs_rule_index, dynmaps);
    if (FALSE == status)
    {
      // Molecule not spoken for yet.  See if it matches.
      status = dynrules_compare_rules(current_lhs_rule_index,
                                      current_rhs_rule_index,
                                      rhs_index,
                                      dynconnections,
                                      dyncolors, dynrules);
      if (TRUE == status)
      {
        rule_t *lhs_rule;
        rule_t *rhs_rule;
        int current_lhs_adjacent_rule_index;
        int current_rhs_adjacent_rule_index;

        // We have a successful match!
        status = dynmaps_map_add(current_lhs_rule_index,
                                 current_rhs_rule_index, dynmaps);
        if (FALSE == status)
        {
          error_printf("Unable to add molecule map.\n");
          return FALSE;
        }

        // Map individual components.
        lhs_rule = dynrules->rules.a + current_lhs_rule_index;
        rhs_rule = dynrules->rules.a + current_rhs_rule_index;
        current_lhs_adjacent_rule_index = lhs_rule->adjacent_index;
        current_rhs_adjacent_rule_index = rhs_rule->adjacent_index;

        while ((-1 != current_lhs_adjacent_rule_index) &&
               (-1 != current_rhs_adjacent_rule_index))
        {
          status = dynmaps_map_add((size_t)current_lhs_adjacent_rule_index,
                                   (size_t)current_rhs_adjacent_rule_index,
                                   dynmaps);
          if (FALSE == status)
          {
            error_printf("Unable to add component map.\n");
            return FALSE;
          }

          lhs_rule = dynrules->rules.a + current_lhs_adjacent_rule_index;
          rhs_rule = dynrules->rules.a + current_rhs_adjacent_rule_index;
          current_lhs_adjacent_rule_index = lhs_rule->and_index;
          current_rhs_adjacent_rule_index = rhs_rule->and_index;
        }

        // Molecule map successfully added.
        return TRUE;
      }
    }

    // Try next rhs molecule.
    current_rhs_rule_index++;
    current_rhs_rule++;
  }

  // Didn't match any of the rhs molecules.
  error_printf("lhs molecule didn't match any rhs molecules.");
  return FALSE;
}

boolean
dynrules_compare_rules(size_t lhs_rule_index, size_t rhs_rule_index,
                       size_t rhs_index,
                       dynconnections_t *dynconnections,
                       dyncolors_t *dyncolors, dynrules_t *dynrules)
{
  boolean status;
  rule_t *lhs_rule;
  rule_t *rhs_rule;
  rule_t *lhs_component_rule;
  rule_t *rhs_component_rule;

  if (lhs_rule_index == -1)
  {
    // lhs rule does not exist.
    if (rhs_rule_index == -1)
    {
      // lhs rule and rhs rule do not exist, so they match.
      return TRUE;
    }
    else
    {
      // lhs rules does not exist, but rhs rule exists, so no match.
      return FALSE;
    }
  }
  else if (rhs_rule_index == -1)
  {
    // lhs rule exists but rhs rule does not exist, so no match.
    return FALSE;
  }


  lhs_rule = dynrules->rules.a + lhs_rule_index;
  rhs_rule = dynrules->rules.a + rhs_rule_index;

  // Check that color is the same.
  if (lhs_rule->color_value != rhs_rule->color_value)
  {
    return FALSE;
  }

  // Check each component in turn.

  /*
   * Rules are being compared so that connections can be applied.
   * Counts need to be * matched up so that appropriate rules are
   * paired together.
   * Verify that each component is appropriately present.
   * States may not be the same (if state is meant to change), so simply test
   * for existance of state.
   */
  // Check that count is the same.
  if ((lhs_rule->count != rhs_rule->count) ||
      (lhs_rule->count_flag != rhs_rule->count_flag))
  {
    return FALSE;
  }


  for (lhs_component_rule = lhs_rule->adjacent_rule,
       rhs_component_rule = rhs_rule->adjacent_rule;
       ;
       lhs_component_rule = lhs_component_rule->and_rule,
       rhs_component_rule = rhs_component_rule->and_rule)
  {
    if (NULL == lhs_component_rule)
    {
      if (NULL == rhs_component_rule)
      {
        // No more components to check.
        break;
      }
      else
      {
        // Component number does not match up.
        return FALSE;
      }
    }
    else if (NULL == rhs_component_rule)
    {
      // Component number does not match up.
      return FALSE;
    }

    // Check that count is the same.
    if ((lhs_component_rule->count != rhs_component_rule->count) ||
        (lhs_component_rule->count_flag != rhs_component_rule->count_flag))
    {
      return FALSE;
    }

    // Check that any connections are the same.
    status = dynconnections_compare_rules(lhs_component_rule->index,
                                          rhs_component_rule->index,
                                          rhs_index, dynconnections);
    if (FALSE == status)
    {
      return FALSE;
    }
  }

  // Rules matched.
  return TRUE;
}

boolean
dynrules_complete(dynrules_t *dynrules)
{
  int i;
  rule_t *current_rule;

  current_rule = dynrules->rules.a;
  for (i = 0; i < dynrules->count; i++)
  {
    if (0 == rule_complete(current_rule, dynrules))
    {
      // Error message already output by rule_complete.
      assert(0);
      return FALSE;
    }

    current_rule++;
  }

  // Found all rules.
  return TRUE;
}

dynrules_t *
dynrules_create(const size_t init_count)
{
  dynrules_t *new_dynrules;
  int status;

  new_dynrules = calloc(1, sizeof(dynrules_t));
  if (NULL == new_dynrules)
  {
    error_printf("Unable to allocate new dynrules.\n");
    assert(0);
    return NULL;
  }

  status = dynarray_create(init_count,
                           &new_dynrules->alloc_count, sizeof(rule_t),
                           &new_dynrules->rules.v);
  if (0 == status)
  {
    // Error message output in dynarray_create.
    free(new_dynrules);
    assert(0);
    return NULL;
  }

  new_dynrules->count = 0;

  return new_dynrules;
}

static boolean
dynrules_resize(dynrules_t *dynrules, const size_t count)
{
  int status;

  status = dynarray_resize(count, &dynrules->alloc_count, sizeof(rule_t),
                           &dynrules->rules.v);
  if (0 == status)
  {
    // Error message output in dynarray_resize.
    return FALSE;
  }

  return TRUE;
}

void
dynrules_output_rule(const output_type ptype, const int level,
                     rule_t *rule, dynrules_t *dynrules,
                     dyncolors_t *dyncolors)
{
  size_t count;
  boolean *used;

  count = dynrules->count;

  used = calloc(sizeof(boolean), count);
  if (NULL == used)
  {
    assert(0);
    return;
  }

  rule_output(ptype, level, rule, dyncolors, "  ", used);

  free(used);
}

void
dynrules_output(const output_type ptype, const int level,
                dynrules_t *dynrules, dyncolors_t *dyncolors)
{
  int i;
  rule_t *current_rule;
  size_t count;
  boolean *used;

  count = dynrules->count;

  used = calloc(sizeof(boolean), count);
  if (NULL == used)
  {
    assert(0);
    return;
  }

  for (i = 0, current_rule = dynrules->rules.a;
       i < count;
       i++, current_rule++)
  {
    rule_output(ptype, level, current_rule, dyncolors, "", used);
  }

  free(used);
}

void
dynrules_destroy(dynrules_t *dynrules)
{
  int i;
  rule_t *current_rule;

  for (i = 0, current_rule = dynrules->rules.a;
       i < dynrules->count;
       i++, current_rule++)
  {
    rule_destroy(current_rule);
  }

  free(dynrules->rules.a);
  free(dynrules);
}
