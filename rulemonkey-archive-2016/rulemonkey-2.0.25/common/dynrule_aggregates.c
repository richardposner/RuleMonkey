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


#include "rule_aggregate.h"
#include "dynrule_aggregates.h"
#include "dynarray.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "constants.h"
#include "output.h"

void suppress_warnings12() {SUPPRESS_WARNINGS;}

static int dynrule_aggregates_get_bind_nodes
             (dynrule_aggregates_t *a, set *set_a);

/**
** Find the number of nodes that are:
** a_onlyp : only bound in first reactant
** b_onlyp : only bound in second reactant
** ab_bothp : bound in 1st reactant and bound in 2nd reactant
*/
int
dynrule_aggregates_count_shared_components
  (dynrule_aggregates_t *a, dynrule_aggregates_t *b,
   specie_t *specie, dyncolors_t *dyncolors,
   double *a_onlyp, double *b_onlyp, double *ab_bothp)
{
  set *set_a;
  set *set_b;
  graph *g;
  int m;
  int a_only;
  int b_only;
  int ab_both;
  int status;
  int i;
  int *colors;

  colors = specie->cg->colors.a;
  g = specie->cg->g.a;
  m = specie->cg->m;

  set_a = calloc(specie->cg->m, sizeof(set));
  set_b = calloc(specie->cg->m, sizeof(set));

  // Get set from a.
  status = dynrule_aggregates_get_bind_nodes(a, set_a);
  if (0 == status)
  {
    free(set_a);
    free(set_b);

    return 0;
  }

  // Get set from b.
  status = dynrule_aggregates_get_bind_nodes(b, set_b);
  if (0 == status)
  {
    free(set_a);
    free(set_b);

    return 0;
  }

  // For each node in a, see if it exists in b.
  a_only = 0;
  b_only = 0;
  ab_both = 0;
  for (i = -1; (i = nextelement(set_a, m, i)) >= 0;)
  {
    if (ISELEMENT(set_b,i))
    {
      // In both.
      ab_both++;
      // Delete element so it isn't counted more than once.
      DELELEMENT(set_b,i);
    }
    else
    {
      // Just in a.
      a_only++;
    }
  }

  // Count nodes that are only in b.
  for (i = -1; (i = nextelement(set_b, m, i)) >= 0;)
  {
    // Node was not found in set a.
    b_only++;
  }

  *a_onlyp = a_only;
  *b_onlyp = b_only;
  *ab_bothp = ab_both;

  free(set_a);
  free(set_b);

  return 1;
}

/**
** For each rule aggregate in a, add element to set_a that is being bound.
*/
static int
dynrule_aggregates_get_bind_nodes(dynrule_aggregates_t *a, set *set_a)
{
  int i;
  int *current_aggregate;
  int status;

  for (i = 0, current_aggregate =
              a->rule_aggregates.a;
       i < a->count;
       i++, current_aggregate += a->rule_aggregate_num_ints)
  {
    int node;

    status = rule_aggregate_first_bind(current_aggregate, &node);
    if (-1 == node)
    {
      return 0;
    }

    // Add found element to set.
    ADDELEMENT(set_a, node);
  }

  return 1;
}

int
dynrule_aggregates_identical
  (dynrule_aggregates_t *a, dynrule_aggregates_t *b)
{
  int status;
  int i;
  size_t agg_size;

  agg_size = a->rule_aggregate_size;

  if (agg_size != b->rule_aggregate_size)
  {
    // No match
    return 0;
  }

  if (a->count != b->count)
  {
    // No match.
    return 0;
  }

  for (i = 0; i < a->count; i++)
  {
    status = rule_aggregate_compare(a->rule_aggregates.a, b->rule_aggregates.a);
    if (status != 0)
    {
      // No match.
      return 0;
    }
  }

  // All matches succeed, so identical.
  return 1;
}

size_t
dynrule_aggregates_sizeof_full_rule_aggregate
  (dynrule_aggregates_t *dynrule_aggregates)
{
  return rule_aggregate_sizeof_full
                (dynrule_aggregates->split_count,
                 dynrule_aggregates->bind_count,
                 dynrule_aggregates->self_bind_count,
                 dynrule_aggregates->state_change_count,
                 dynrule_aggregates->transient_bind_count,
                 dynrule_aggregates->transient_self_bind_count);
}

void
dynrule_aggregates_output(const output_type ptype, const int level,
                          dynrule_aggregates_t *dynrule_aggregates)
{
  size_t i;
  int *current_aggregate;

  custom_printf(ptype, level,
         "dynrule_aggregates: s:%zu b:%zu sb:%zu sc:%zu tb:%zu tsb:%zu; "
         "num_ints:%zu count:%zu\n",
         dynrule_aggregates->split_count,
         dynrule_aggregates->bind_count, dynrule_aggregates->self_bind_count,
         dynrule_aggregates->state_change_count,
         dynrule_aggregates->transient_bind_count,
         dynrule_aggregates->transient_self_bind_count,
         dynrule_aggregates->rule_aggregate_num_ints,
         dynrule_aggregates->count);
  custom_printf(ptype, level, "All aggregates:\n");
  for (i = 0, current_aggregate =
              dynrule_aggregates->rule_aggregates.a;
       i < dynrule_aggregates->count;
       i++, current_aggregate += dynrule_aggregates->rule_aggregate_num_ints)
  {
    rule_aggregate_output(ptype, level, current_aggregate);
    custom_printf(ptype, level, "\n");
  }
}

dynrule_aggregates_t *
dynrule_aggregates_create(const size_t init_count,
                          const size_t split_count,
                          const size_t bind_count,
                          const size_t self_bind_count,
                          const size_t state_change_count,
                          const size_t transient_bind_count,
                          const size_t transient_self_bind_count,
                          const int product_count)
{
  int status;
  int *current;

  dynrule_aggregates_t *new_dynrule_aggregates;
  size_t current_init_count;

  // Allocate new rule aggregate.
  new_dynrule_aggregates = calloc((size_t)1, sizeof(dynrule_aggregates_t));
  if (NULL == new_dynrule_aggregates)
  {
    error_printf("Allocating new rule aggregate");
    assert(0);
    return NULL;
  }

  // Set structure members.
  new_dynrule_aggregates->split_count = split_count;
  new_dynrule_aggregates->bind_count = bind_count;
  new_dynrule_aggregates->self_bind_count = self_bind_count;
  new_dynrule_aggregates->state_change_count = state_change_count;
  new_dynrule_aggregates->transient_bind_count = transient_bind_count;
  new_dynrule_aggregates->transient_self_bind_count = transient_self_bind_count;
  new_dynrule_aggregates->product_count = product_count;

  new_dynrule_aggregates->rule_aggregate_num_ints =
    rule_aggregate_num_ints (split_count, bind_count, self_bind_count,
                             state_change_count, transient_bind_count,
                             transient_self_bind_count);
  new_dynrule_aggregates->rule_aggregate_size =
    new_dynrule_aggregates->rule_aggregate_num_ints * sizeof(int);

  // Always create at least one so that some memory is allocated.
  if (init_count < 1)
  {
    current_init_count = 1;
  }
  else
  {
    current_init_count = init_count;
  }

  status = dynarray_create
             (current_init_count,
              &new_dynrule_aggregates->alloc_count,
              new_dynrule_aggregates->rule_aggregate_size,
              &new_dynrule_aggregates->rule_aggregates.v);
  if (0 == status)
  {
    // Error message output in dynarrays_create.
    free(new_dynrule_aggregates);
    return NULL;
  }

  // Allocate memory to hold template
  new_dynrule_aggregates->rule_aggregate_template =
    calloc(1, new_dynrule_aggregates->rule_aggregate_size);

  // Initialize terminators to delimit block types.
  current = new_dynrule_aggregates->rule_aggregate_template;
  *current++ = 1;
  current += split_count * AGG_SPLIT_LENGTH;
  *current++ = -1;
  current += bind_count * AGG_BIND_LENGTH;
  *current++ = -1;
  current += self_bind_count * AGG_SELF_BIND_LENGTH;
  *current++ = -1;
  current += state_change_count * AGG_STATE_CHANGE_LENGTH;
  *current++ = -1;
  current += transient_bind_count * AGG_BIND_LENGTH;
  *current++ = -1;
  current += transient_self_bind_count * AGG_SELF_BIND_LENGTH;
  *current = -1;

  // Nothing allocated yet.
  new_dynrule_aggregates->count = 0;

  return new_dynrule_aggregates;
}

int
dynrule_aggregates_rule_aggregate_create
  (dynrule_aggregates_t *dynrule_aggregates)
{
  int status;
  int *new_aggregate;

  dynrule_aggregates->count++;

  status =
    dynarray_resize(dynrule_aggregates->count,
                    &dynrule_aggregates->alloc_count,
                    dynrule_aggregates->rule_aggregate_size,
                    &dynrule_aggregates->rule_aggregates.v);
  if (0 == status)
  {
    error_printf("dynrule_aggregates_rule_aggregate_create: Resize error.\n");
    return FALSE;
  }

  // Get pointer to memory just allocated.
  new_aggregate =
    dynrule_aggregates->rule_aggregates.a +
      ((dynrule_aggregates->count - 1) *
       dynrule_aggregates->rule_aggregate_num_ints);
  
  // Copy template to initialize new aggregate.
  memcpy(new_aggregate, dynrule_aggregates->rule_aggregate_template,
         dynrule_aggregates->rule_aggregate_size);

  // Save point to new aggregate so it is ready to use.
  dynrule_aggregates->current_aggregate = new_aggregate;

  return TRUE;
}

int
dynrule_aggregates_add_rule(dynrule_aggregates_t *dynrule_aggregates,
                            rule_t *rule, const int matched_node)
{

  int *rule_aggregate;
  int conflict;

  conflict = FALSE;

  rule_aggregate = dynrule_aggregates->current_aggregate;

  if (FALSE == rule->attributes_modified)
  {
    // Don't bother adding rule because rule doesn't contain any modifiers.
#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "No modifiers in rule %d\n", rule->index);
#endif

    return TRUE;
  }

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(8, "\nAdding rule %d.\n", rule->index);
  DYNRULE_AGGREGATES_OUTPUT(OUTPUT_DEBUG, 8, dynrule_aggregates);
#endif

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(3, "Adding rule %d (node %d) ", rule->index, matched_node);
  RULE_AGGREGATE_OUTPUT(OUTPUT_DEBUG, 3, rule_aggregate);
  DEBUG_PRINTF(3, " -> ");
#endif

  // Skip multiplier.
  rule_aggregate++;

  // Split directives.
  rule_aggregate = rule_aggregate_double_add_rule
                (rule_aggregate, dynrule_aggregates->split_count,
                 rule->primary_split_index_list,
                 rule->secondary_split_index_list, matched_node, &conflict);
  if (TRUE == conflict) return FALSE;

  // Bind directives.
  rule_aggregate = rule_aggregate_single_add_rule
                (rule_aggregate, dynrule_aggregates->bind_count,
                 rule->bind_index_list, matched_node, &conflict);
  if (TRUE == conflict) return FALSE;

  // Self bind directives.
  rule_aggregate = rule_aggregate_double_add_rule
                (rule_aggregate, dynrule_aggregates->self_bind_count,
                 rule->primary_self_bind_index_list,
                 rule->secondary_self_bind_index_list, matched_node, &conflict);
  if (TRUE == conflict) return FALSE;

  // State change directives.
  rule_aggregate = rule_aggregate_state_change_add_rule
                (rule_aggregate, dynrule_aggregates->state_change_count,
                 rule->state_change_index_list,
                 rule->state_change_old_list,
                 rule->state_change_new_list,
                 matched_node, &conflict);
  if (TRUE == conflict) return FALSE;

  // Transient bind directives.
  rule_aggregate = rule_aggregate_single_add_rule
                (rule_aggregate, dynrule_aggregates->transient_bind_count,
                 rule->transient_bind_index_list, matched_node, &conflict);
  if (TRUE == conflict) return FALSE;

  // Transient self bind directives.
  rule_aggregate = rule_aggregate_double_add_rule
                (rule_aggregate,
                 dynrule_aggregates->transient_self_bind_count,
                 rule->primary_transient_self_bind_index_list,
                 rule->secondary_transient_self_bind_index_list,
                 matched_node, &conflict);
  if (TRUE == conflict) return FALSE;

#ifdef RULE_CHECK_DEBUG
  RULE_AGGREGATE_OUTPUT(OUTPUT_DEBUG, 3,
                        dynrule_aggregates->current_aggregate);
  DEBUG_PRINTF(3, "\n");
#endif

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(8, "Added rule.\n");
  DYNRULE_AGGREGATES_OUTPUT(OUTPUT_DEBUG, 8, dynrule_aggregates);
#endif

  return TRUE;
}

void
dynrule_aggregates_destroy(dynrule_aggregates_t *dynrule_aggregates)
{
  free(dynrule_aggregates->rule_aggregates.a);
  free(dynrule_aggregates->rule_aggregate_template);

  free(dynrule_aggregates);
}

void
dynrule_aggregates_reset_all(dynrule_aggregates_t *dynrule_aggregates)
{
  if (NULL == dynrule_aggregates)
  {
    // Nothing to do.
    return;
  }

  dynrule_aggregates->count = 0;
}

/**
 * Check if the current aggregate is a duplicate of a previously
 * added aggregate.  If so, remove aggregate from list.
 */
int
dynrule_aggregates_check_current_unique
  (dynrule_aggregates_t *dynrule_aggregates)
{
  int i;
  int count;
  int *existing_aggregate;
  int *current_aggregate;

  // Don't check last aggregate against itself;
  count = dynrule_aggregates->count - 1;

  current_aggregate = dynrule_aggregates->current_aggregate;

  for (i = 0, existing_aggregate =
              dynrule_aggregates->rule_aggregates.a;
       i < count;
       i++, existing_aggregate += dynrule_aggregates->rule_aggregate_num_ints)
  {
    if (1 == rule_aggregate_compare(existing_aggregate, current_aggregate))
    {
#ifdef RULE_CHECK_DEBUG
      DEBUG_PRINTF(4, "Aggregate already exists: ");
      RULE_AGGREGATE_OUTPUT(OUTPUT_DEBUG, 4, current_aggregate);
      DEBUG_PRINTF(4, "\n");
#endif

      // Found a match, so remove last aggregate
      dynrule_aggregates->count = count;
      dynrule_aggregates->current_aggregate -=
        dynrule_aggregates->rule_aggregate_num_ints;

      return 0;
    }
  }

  // No match was found.
  return 1;
}

int
dynrule_aggregates_check_product_count
      (dynrule_aggregates_t *dynrule_aggregates, specie_t *specie)
{
  int product_count;

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(8, "\nChecking product count\n");
  DYNRULE_AGGREGATES_OUTPUT(OUTPUT_DEBUG, 8, dynrule_aggregates);
#endif

  // Calculate number of products produced by current aggregate
  product_count = rule_aggregate_calculate_product_count
            (specie->cg, dynrule_aggregates->current_aggregate);

  // Compare current product count with correct product count.
  if (dynrule_aggregates->product_count != product_count)
  {
#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(4, "Product count: %d does not match expected: %d.\n",
                 product_count, dynrule_aggregates->product_count);
    RULE_AGGREGATE_OUTPUT(OUTPUT_DEBUG, 4,
                          dynrule_aggregates->current_aggregate);
    DEBUG_PRINTF(4, "\n");
#endif

    // Product count not correct, so remove last aggregate.
    dynrule_aggregates->count--;
    dynrule_aggregates->current_aggregate -=
      dynrule_aggregates->rule_aggregate_num_ints;

    // Current product count does not equal correct product count,
    // so reaction rule does not match.
    return 0;
  }

  return 1;
}

