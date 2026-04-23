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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "constants.h"
#include "output.h"
#include "cgraph.h"

void suppress_warnings38() {SUPPRESS_WARNINGS;}

static void
rule_aggregate_raw_add_rule(int *rule_aggregate,
                            const size_t used_offset,
                            const size_t offset,
                            const size_t sec_used_offset,
                            const size_t block_length,
                            const ilist_t *current_node,
                            const int matched_node, int *changed);

static boolean
rule_aggregate_double_first_active(int *ppri, int *psec, int **prule_aggregate);

static boolean
rule_aggregate_single_first_active(int *ppri, int *psec, int **prule_aggregate);

static boolean
rule_aggregate_state_change_first_active(int *ppri, int *psec,
                                         int **prule_aggregate);

static boolean
rule_aggregate_state_change_merge_partial(int *base_full, int *base_b,
                                          void *data_full, void *data_b);

static boolean
rule_aggregate_raw_merge_partial(int *base_full, int *base_b,
                                 int *renum_full, int *renum_b,
                                 const int full_used_offset,
                                 const int full_node_offset,
                                 const int b_used_offset,
                                 const int b_node_offset);

static boolean
rule_aggregate_single_merge_partial(int *base_full, int *base_b,
                                    void *renum_full, void *renum_b);

static boolean
rule_aggregate_double_merge_partial(int *base_full, int *base_b,
                                    void *renum_full, void *renum_b);

static boolean
rule_aggregate_single_cgraph_append(int **pbase_full, int **pbase_b,
                                    int *renum_a, int *renum_b);

static boolean
rule_aggregate_state_change_cgraph_append(int **pbase_full, int **pbase_b,
                                          int *renum_a, int *renum_b);

static int
rule_aggregate_split(int *base_a, void *data1);

static int
rule_aggregate_bind(int *base_a, void *data1);

static int
rule_aggregate_state_change(int *base_a, void *data1);

typedef int (*rule_aggregate_iterate_cb_t)(int *base_a, void *data1);
typedef boolean (*rule_aggregate_iterate2_cb_t)(int *base_a, int *base_b,
                                                void *data1, void *data2);
static boolean
rule_aggregate_skip(int **pbase_a, const size_t element_a_size);

static boolean
rule_aggregate_iterate(rule_aggregate_iterate_cb_t rule_aggregate_iterate_cb,
                       int **pbase_a,
                       void *data1,
                       const size_t element_a_size);

static boolean
rule_aggregate_iterate2(rule_aggregate_iterate2_cb_t rule_aggregate_iterate2_cb,
                       int **pbase_a, int **pbase_b,
                       const size_t element_a_length,
                       const size_t element_b_length,
                       void *data1, void *data2);

static int *
rule_aggregate_output_single(const output_type ptype, const int level,
                             int *base);

static int *
rule_aggregate_output_double(const output_type ptype, const int level,
                             int *base);

static int *
rule_aggregate_output_state_change(const output_type ptype, const int level,
                                   int *base);





static boolean
rule_aggregate_skip(int **pbase_a, const size_t element_a_size)
{
  int *base_a;

  base_a = *pbase_a;

  // Iterate through blocks.
  while(*base_a != -1)
  {
    base_a += element_a_size;
  }

  // Skip terminator.
  base_a++;

  *pbase_a = base_a;

  return TRUE;
}

static boolean
rule_aggregate_iterate(rule_aggregate_iterate_cb_t rule_aggregate_iterate_cb,
                       int **pbase_a,
                       void *data1,
                       const size_t element_a_size)
{
  int *base_a;
  int status;

  base_a = *pbase_a;

  // Iterate through blocks.
  while(*base_a != -1)
  {
    status = rule_aggregate_iterate_cb(base_a, data1);
    if (FALSE == status) return FALSE;

    base_a += element_a_size;
  }

  // Skip terminator.
  base_a++;

  *pbase_a = base_a;

  return TRUE;
}

static boolean
rule_aggregate_iterate2(rule_aggregate_iterate2_cb_t rule_aggregate_iterate2_cb,
                       int **pbase_a, int **pbase_b,
                       const size_t element_a_length,
                       const size_t element_b_length,
                       void *data1, void *data2)
{
  int *base_a;
  int *base_b;
  boolean status;

  base_a = *pbase_a;
  base_b = *pbase_b;

  // Iterate through blocks.
  while(*base_a != -1)
  {
    assert(-1 != *base_b);

    status = rule_aggregate_iterate2_cb(base_a, base_b, data1, data2);
    if (FALSE == status) return FALSE;

    base_a += element_a_length;
    base_b += element_b_length;
  }

  // Skip terminator.
  assert(-1 == *base_b);
  base_a++;
  base_b++;

  *pbase_a = base_a;
  *pbase_b = base_b;

  return TRUE;
}

static boolean
rule_aggregate_raw_aggregate_compare
  (const int **pbase1, const int **pbase2,
   const size_t used_offset, const size_t node_offset,
   const size_t block_length);

static boolean
rule_aggregate_double_aggregate_compare
  (const int **pbase1, const int **pbase2);

static boolean
rule_aggregate_state_change_aggregate_compare
  (const int **pbase1, const int **pbase2);

static boolean
rule_aggregate_single_aggregate_compare
  (const int **pbase1, const int **pbase2);



size_t
rule_aggregate_sizeof_full(const size_t split_count,
                           const size_t bind_count,
                           const size_t self_bind_count,
                           const size_t state_change_count,
                           const size_t transient_bind_count,
                           const size_t transient_self_bind_count)
{
  return (1 +
          (split_count * FULL_AGG_SPLIT_LENGTH) + 1 +
          (bind_count * FULL_AGG_BIND_LENGTH) + 1 +
          (self_bind_count * FULL_AGG_SELF_BIND_LENGTH) + 1 +
          (state_change_count * FULL_AGG_STATE_CHANGE_LENGTH) + 1 +
          (transient_bind_count * FULL_AGG_BIND_LENGTH) + 1 +
          (transient_self_bind_count * FULL_AGG_SELF_BIND_LENGTH) + 1
         ) *
         sizeof(int);
}

size_t
rule_aggregate_num_ints(const size_t split_count,
                        const size_t bind_count,
                        const size_t self_bind_count,
                        const size_t state_change_count,
                        const size_t transient_bind_count,
                        const size_t transient_self_bind_count)
{
  return (1 +
          (split_count * AGG_SPLIT_LENGTH) + 1 +
          (bind_count * AGG_BIND_LENGTH) + 1 +
          (self_bind_count * AGG_SELF_BIND_LENGTH) + 1 +
          (state_change_count * AGG_STATE_CHANGE_LENGTH) + 1 +
          (transient_bind_count * AGG_BIND_LENGTH) + 1 +
          (transient_self_bind_count * AGG_SELF_BIND_LENGTH) + 1
         );
}

void
rule_aggregate_init_full(int *full_rule_aggregate, int *rule_aggregate)
{
  // Initialize multiplier.
  full_rule_aggregate[AGG_MULTIPLIER_OFFSET] = 1;

  // Skip multiplier.
  rule_aggregate++;
  full_rule_aggregate++;

  // Split.
  while (*rule_aggregate != -1)
  {
    memset(full_rule_aggregate, 0, FULL_AGG_SPLIT_LENGTH * sizeof(int));
    rule_aggregate += AGG_SPLIT_LENGTH;
    full_rule_aggregate += FULL_AGG_SPLIT_LENGTH;
  }
  // Add delimiter.
  *full_rule_aggregate = -1;
  // Skip delimiter.
  rule_aggregate++;
  full_rule_aggregate++;

  // Bind.
  while (*rule_aggregate != -1)
  {
    memset(full_rule_aggregate, 0, FULL_AGG_BIND_LENGTH * sizeof(int));
    rule_aggregate += AGG_BIND_LENGTH;
    full_rule_aggregate += FULL_AGG_BIND_LENGTH;
  }
  // Add delimiter.
  *full_rule_aggregate = -1;
  // Skip delimiter.
  rule_aggregate++;
  full_rule_aggregate++;

  // Self bind.
  while (*rule_aggregate != -1)
  {
    memset(full_rule_aggregate, 0, FULL_AGG_SELF_BIND_LENGTH * sizeof(int));
    rule_aggregate += AGG_SELF_BIND_LENGTH;
    full_rule_aggregate += FULL_AGG_SELF_BIND_LENGTH;
  }
  // Add delimiter.
  *full_rule_aggregate = -1;
  // Skip delimiter.
  rule_aggregate++;
  full_rule_aggregate++;

  // State change.
  while (*rule_aggregate != -1)
  {
    memset(full_rule_aggregate, 0, FULL_AGG_STATE_CHANGE_LENGTH * sizeof(int));
    rule_aggregate += AGG_STATE_CHANGE_LENGTH;
    full_rule_aggregate += FULL_AGG_STATE_CHANGE_LENGTH;
  }
  // Add delimiter.
  *full_rule_aggregate = -1;
  // Skip delimiter.
  rule_aggregate++;
  full_rule_aggregate++;

  // Transient bind.
  while (*rule_aggregate != -1)
  {
    memset(full_rule_aggregate, 0, FULL_AGG_BIND_LENGTH * sizeof(int));
    rule_aggregate += AGG_BIND_LENGTH;
    full_rule_aggregate += FULL_AGG_BIND_LENGTH;
  }
  // Add delimiter.
  *full_rule_aggregate = -1;
  // Skip delimiter.
  rule_aggregate++;
  full_rule_aggregate++;

  // Transient self bind.
  while (*rule_aggregate != -1)
  {
    memset(full_rule_aggregate, 0, FULL_AGG_SELF_BIND_LENGTH * sizeof(int));
    rule_aggregate += AGG_SELF_BIND_LENGTH;
    full_rule_aggregate += FULL_AGG_SELF_BIND_LENGTH;
  }
  // Add delimiter.
  *full_rule_aggregate = -1;
}

static void
rule_aggregate_raw_add_rule(int *rule_aggregate,
                            const size_t used_offset,
                            const size_t offset,
                            const size_t sec_used_offset,
                            const size_t block_length,
                            const ilist_t *current_node,
                            const int matched_node, int *pconflict)
{
  for (;
       NULL != current_node;
       current_node = current_node->next)
  {
    int *base;

    // Get base.
    base = rule_aggregate + current_node->index * block_length;

    if (TRUE == base[used_offset])
    {
      // Element already initialized by previous rule.
      if (FALSE == *pconflict)
      {
        // Report error.
        assert(0);
        *pconflict = TRUE;
        return;
      }
    }
    else
    {
      // Element not used yet.
      // Set to used.
      base[used_offset] = TRUE;
    }

    // Set matched node.
    base[offset] = matched_node;
  }
}

int *
rule_aggregate_state_change_add_rule(int *rule_aggregate,
                                     const size_t count,
                                     ilist_t *current_node,
                                     ilist_t *old_list,
                                     ilist_t *new_list,
                                     const int matched_node, int *pconflict)
{
  for (;
       NULL != current_node && NULL != old_list && NULL != new_list;
       current_node = current_node->next,
       old_list = old_list->next, new_list = new_list->next)
  {
    int *base;

    // Get base.
    base = rule_aggregate + current_node->index * AGG_STATE_CHANGE_LENGTH;

    if (TRUE == base[AGG_USED_OFFSET])
    {
      // Element already initialized by previous rule.
      if (FALSE == *pconflict)
      {
        // Report error.
        assert(0);
        *pconflict = TRUE;
        return NULL;
      }
    }
    else
    {
      // Element not used yet.
      // Set to used.
      base[AGG_USED_OFFSET] = TRUE;
    }

    // Set matched node.
    base[AGG_NODE_OFFSET] = matched_node;
    base[AGG_OLD_STATE_OFFSET] = old_list->index;
    base[AGG_NEW_STATE_OFFSET] = new_list->index;
  }

  // Skip all directives and terminator.
  rule_aggregate += count * AGG_STATE_CHANGE_LENGTH + 1;

#ifndef DS_NO_DOUBLECHECK
  // Already skipped terminator so look back.
  assert(-1 == rule_aggregate[-1]);
#endif

  return  rule_aggregate;
}

int *
rule_aggregate_single_add_rule(int *rule_aggregate,
                               const size_t count,
                               ilist_t *current_node,
                               const int matched_node, int *pconflict)
{
  rule_aggregate_raw_add_rule(rule_aggregate,
                              AGG_PRI_USED_OFFSET,
                              AGG_PRI_NODE_OFFSET,
                              0,
                              AGG_SINGLE_LENGTH,
                              current_node, matched_node, pconflict);

  // Skip all directives and terminator.
  rule_aggregate += count * AGG_SINGLE_LENGTH + 1;

#ifndef DS_NO_DOUBLECHECK
  // Already skipped terminator so look back.
  assert(-1 == rule_aggregate[-1]);
#endif

  return  rule_aggregate;
}

int *
rule_aggregate_double_add_rule(int *rule_aggregate,
                               const size_t count,
                               const ilist_t *current_pri_node,
                               const ilist_t *current_sec_node,
                               const int matched_node, int *pconflict)
{
  int old_conflict;

  old_conflict = *pconflict;

  // Primary directives.
  rule_aggregate_raw_add_rule(rule_aggregate,
                              AGG_PRI_USED_OFFSET,
                              AGG_PRI_NODE_OFFSET,
                              AGG_SEC_USED_OFFSET,
                              AGG_DOUBLE_LENGTH,
                              current_pri_node, matched_node, pconflict);

  if (old_conflict == *pconflict)
  {
    // Conflict not found with previous call, so we can continue.

    // Secondary directives.
    rule_aggregate_raw_add_rule(rule_aggregate,
                                AGG_SEC_USED_OFFSET,
                                AGG_SEC_NODE_OFFSET,
                                0,
                                AGG_DOUBLE_LENGTH,
                                current_sec_node, matched_node, pconflict);
  }

  // Skip all directives and terminator.
  rule_aggregate += count * AGG_DOUBLE_LENGTH + 1;

#ifndef DS_NO_DOUBLECHECK
  // Already skipped terminator so look back.
  assert(-1 == rule_aggregate[-1]);
#endif

  return rule_aggregate;
}

/* If TRUE is returned, contents of prule_aggregate undefined. */
static boolean
rule_aggregate_raw_first_active(const size_t length,
                                int *ppri, int *psec, int **prule_aggregate)
{
  int *rule_aggregate;
  boolean status;

  rule_aggregate = *prule_aggregate;

  status = FALSE;

  // Return first used point.
  for(;;)
  {
    if (-1 == *rule_aggregate)
    {
      // At end of block and no active node found.  Skip terminator.
      rule_aggregate++;
      status = FALSE;
      break;
    }
    else if (FALSE != rule_aggregate[AGG_PRI_USED_OFFSET])
    {
      *ppri = rule_aggregate[AGG_PRI_NODE_OFFSET];
      *psec = -1;

      status = TRUE;
      break;
    }

    rule_aggregate += length;
  }

  if (TRUE == status)
  {
    // Skip over remaining blocks.
    while (*rule_aggregate != -1)
    {
      rule_aggregate += length;
    }
  }

  *prule_aggregate = rule_aggregate;

  return status;
}

/* If TRUE is returned, contents of prule_aggregate undefined. */
static boolean
rule_aggregate_double_first_active(int *ppri, int *psec, int **prule_aggregate)
{
  int *rule_aggregate;
  int status;

  rule_aggregate = *prule_aggregate;

  status = FALSE;

  // Return first used point.
  for(;;)
  {
    if (-1 == *rule_aggregate)
    {
      // At end of block and no active node found.  Skip terminator.
      rule_aggregate++;
      status = FALSE;
      break;
    }
    else if (FALSE != rule_aggregate[AGG_PRI_USED_OFFSET])
    {
      // Return first pair of points.
#ifndef DS_NO_DOUBLECHECK
      assert(FALSE != rule_aggregate[AGG_SEC_USED_OFFSET]);
#endif
      *ppri = rule_aggregate[AGG_PRI_NODE_OFFSET];
      *psec = rule_aggregate[AGG_SEC_NODE_OFFSET];

      status = TRUE;
      break;
    }

    rule_aggregate += AGG_DOUBLE_LENGTH;
  }

  if (TRUE == status)
  {
    // Skip over remaining blocks.
    while (*rule_aggregate != -1)
    {
      rule_aggregate += AGG_DOUBLE_LENGTH;
    }

    // Skip terminator.
    rule_aggregate++;
  }

  *prule_aggregate = rule_aggregate;

  return status;
}

/* If TRUE is returned, contents of prule_aggregate undefined. */
static boolean
rule_aggregate_single_first_active(int *ppri, int *psec, int **prule_aggregate)
{
  return rule_aggregate_raw_first_active(AGG_SINGLE_LENGTH,
                                         ppri, psec, prule_aggregate);
}

/* If TRUE is returned, contents of prule_aggregate undefined. */
static boolean
rule_aggregate_state_change_first_active(int *ppri, int *psec,
                                         int **prule_aggregate)
{
  return rule_aggregate_raw_first_active(AGG_STATE_CHANGE_LENGTH,
                                         ppri, psec, prule_aggregate);
}

boolean
rule_aggregate_first_active(int *aggregate, int *ppri, int *psec,
                            int *pid, boolean *pbimolecular)
{
  int status;

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(6, "Checking first active:\n");
  rule_aggregate_output(OUTPUT_DEBUG, 8, aggregate);
  DEBUG_PRINTF(6, "\n");
#endif

  // Skip multiplier.
  aggregate++;

  // Split directives.
  status = rule_aggregate_double_first_active(ppri, psec, &aggregate);
  if (TRUE == status)
  {
    *pid = *ppri;
    *pbimolecular = FALSE;

    return TRUE;
  }

  // Bind directives.
  status = rule_aggregate_single_first_active(ppri, psec, &aggregate);
  if (TRUE == status)
  {
    *pid = -1;
    *pbimolecular = TRUE;

    return TRUE;
  }

  // Self bind directives.
  status = rule_aggregate_double_first_active(ppri, psec, &aggregate);
  if (TRUE == status)
  {
    *pid = *ppri;
    *pbimolecular = TRUE;

    return TRUE;
  }

  // Status change directives.
  status = rule_aggregate_state_change_first_active(ppri, psec, &aggregate);
  if (TRUE == status)
  {
    *pid = -1;
    *pbimolecular = FALSE;

    return TRUE;
  }

  // Transient bind directives.
  status = rule_aggregate_single_first_active(ppri, psec, &aggregate);
  if (TRUE == status)
  {
    *pid = -1;
    *pbimolecular = TRUE;

    return TRUE;
  }

  // Transient self bind directives.
  status = rule_aggregate_double_first_active(ppri, psec, &aggregate);
  if (TRUE == status)
  {
    *pid = *ppri;
    *pbimolecular = TRUE;

    return TRUE;
  }

  return FALSE;
}

boolean
rule_aggregate_first_bind(int *aggregate, int *ppri)
{
  int status;
  int sec;

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(6, "Checking first active:\n");
  rule_aggregate_output(OUTPUT_DEBUG, 8, aggregate);
  DEBUG_PRINTF(6, "\n");
#endif

  // Skip multiplier.
  aggregate++;

  // Split directives.
  rule_aggregate_skip(&aggregate, AGG_DOUBLE_LENGTH);

  // Bind directives.
  status = rule_aggregate_single_first_active(ppri, &sec, &aggregate);
  if (TRUE == status)
  {
    return TRUE;
  }

  // Skip self bind directives.
  rule_aggregate_skip(&aggregate, AGG_DOUBLE_LENGTH);

  // Skip status change directives.
  rule_aggregate_skip(&aggregate, AGG_STATE_CHANGE_LENGTH);

  // Transient bind directives.
  status = rule_aggregate_single_first_active(ppri, &sec, &aggregate);
  if (TRUE == status)
  {
    return TRUE;
  }

  error_printf("Unable to find first bind (self bind doesn't count).\n");

  return FALSE;
}

/*
** Merge state change partial rule_aggregate block into a full rule_aggregate block.
*/
static boolean
rule_aggregate_state_change_merge_partial(int *base_full, int *base_b,
                                          void *data_full, void *data_b)
{
  int *renum_full;
  int *renum_b;

  renum_full = data_full;
  renum_b = data_b;

  if (TRUE == base_b[AGG_PRI_USED_OFFSET])
  {
    if (FALSE == base_full[AGG_PRI_USED_OFFSET])
    {
      base_full[AGG_PRI_USED_OFFSET] = TRUE;

      if (NULL == renum_b)
      {
        // Add b.
        base_full[AGG_PRI_NODE_OFFSET] = base_b[AGG_PRI_NODE_OFFSET];
      }
      else
      {
        // Add b, renumbering node appropriately.
        base_full[AGG_PRI_NODE_OFFSET] = renum_b[base_b[AGG_PRI_NODE_OFFSET]];
      }
      base_full[AGG_OLD_STATE_OFFSET] = base_b[AGG_OLD_STATE_OFFSET];
      base_full[AGG_NEW_STATE_OFFSET] = base_b[AGG_NEW_STATE_OFFSET];
    }
    else
    {
      // point already spoken for.
      assert(0);
      return FALSE;
    }
  }
  else if (TRUE == base_full[AGG_PRI_USED_OFFSET])
  {
    // Need to renumber full rule_aggregate.
    base_full[AGG_PRI_NODE_OFFSET] = renum_full[base_full[AGG_PRI_NODE_OFFSET]];
  }

  return TRUE;
}

static boolean
rule_aggregate_raw_merge_partial(int *base_full, int *base_b,
                                 int *renum_full, int *renum_b,
                                 const int full_used_offset,
                                 const int full_node_offset,
                                 const int b_used_offset,
                                 const int b_node_offset)
{
  if (TRUE == base_b[b_used_offset])
  {
    // Something from b to add.
    if (FALSE == base_full[full_used_offset])
    {
      // Point not used yet.
      base_full[full_used_offset] = TRUE;

      if (NULL == renum_b)
      {
        // Add b.
        base_full[full_node_offset] = base_b[b_node_offset];
      }
      else
      {
        // Add b, renumbering node appropriately.
        base_full[full_node_offset] = renum_b[base_b[b_node_offset]];
      }
    }
    else
    {
      // point already spoken for.
      // Need to renumber full rule_aggregate.
      assert(NULL != renum_full);
      base_full[full_node_offset] = renum_full[base_full[full_node_offset]];
      return FALSE;
    }
  }
  else if (TRUE == base_full[full_used_offset])
  {
    // Nothing from b to add.
    // Need to renumber full rule_aggregate.
    base_full[full_node_offset] = renum_full[base_full[full_node_offset]];
  }

  return TRUE;
}

static boolean
rule_aggregate_single_merge_partial(int *base_full, int *base_b,
                                    void *renum_full, void *renum_b)
{
  boolean status;

  status = rule_aggregate_raw_merge_partial(base_full, base_b,
                                            renum_full, renum_b,
                                            AGG_PRI_USED_OFFSET,
                                            AGG_PRI_NODE_OFFSET,
                                            AGG_PRI_USED_OFFSET,
                                            AGG_PRI_NODE_OFFSET);
  if (FALSE == status)
  {
    // First spot is already taken, try adding to second spot.
    status = rule_aggregate_raw_merge_partial(base_full, base_b,
                                              renum_full, renum_b,
                                              AGG_SEC_USED_OFFSET,
                                              AGG_SEC_NODE_OFFSET,
                                              AGG_PRI_USED_OFFSET,
                                              AGG_PRI_NODE_OFFSET);
  }

  return status;
}

static boolean
rule_aggregate_double_merge_partial(int *base_full, int *base_b,
                                    void *renum_full, void *renum_b)
{
  boolean status;

  status = rule_aggregate_raw_merge_partial(base_full, base_b,
                                            renum_full, renum_b,
                                            AGG_PRI_USED_OFFSET,
                                            AGG_PRI_NODE_OFFSET,
                                            AGG_PRI_USED_OFFSET,
                                            AGG_PRI_NODE_OFFSET);
  if (FALSE == status)
  {
    return FALSE;
  }

  status = rule_aggregate_raw_merge_partial(base_full, base_b,
                                            renum_full, renum_b,
                                            AGG_SEC_USED_OFFSET,
                                            AGG_SEC_NODE_OFFSET,
                                            AGG_SEC_USED_OFFSET,
                                            AGG_SEC_NODE_OFFSET);
  
  return status;
}

static boolean
rule_aggregate_single_cgraph_append(int **pbase_full, int **pbase_b,
                                    int *renum_a, int *renum_b)
{
  return rule_aggregate_iterate2
           (rule_aggregate_single_merge_partial, pbase_full, pbase_b,
            FULL_AGG_DOUBLE_LENGTH, AGG_SINGLE_LENGTH, renum_a, renum_b);
}

static boolean
rule_aggregate_double_cgraph_append(int **pbase_full, int **pbase_b,
                                    int *renum_a, int *renum_b)
{
  return rule_aggregate_iterate2
           (rule_aggregate_double_merge_partial, pbase_full, pbase_b,
            FULL_AGG_DOUBLE_LENGTH, AGG_DOUBLE_LENGTH, renum_a, renum_b);
}

static boolean
rule_aggregate_state_change_cgraph_append(int **pbase_full, int **pbase_b,
                                          int *renum_a, int *renum_b)
{
  return rule_aggregate_iterate2
           (rule_aggregate_state_change_merge_partial, pbase_full, pbase_b,
            FULL_AGG_STATE_CHANGE_LENGTH, AGG_STATE_CHANGE_LENGTH,
            renum_a, renum_b);
}

boolean
rule_aggregate_cgraph_append(cgraph_t *a, int *rule_aggregate_full,
                             cgraph_t *b, int *rule_aggregate_b,
                             cgraph_t *c)
{
  int *renum_a; // Maps element numbers from a to c
  int *renum_b; // Maps element numbers from b to c
  int status;

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(5, "In cgraph append\n");
  DEBUG_PRINTF(5, "Aggregate full: ");
  RULE_AGGREGATE_OUTPUT_FULL(OUTPUT_DEBUG, 5, rule_aggregate_full);
  DEBUG_PRINTF(5, "\nAggregate b: ");
  RULE_AGGREGATE_OUTPUT(OUTPUT_DEBUG, 5, rule_aggregate_b);
  DEBUG_PRINTF(5, "\n");
#endif

  if (NULL != a)
  {
    renum_a = calloc((size_t)a->n, sizeof(int));
    if (NULL == renum_a)
    {
      assert(0);
      return FALSE;
    }

    renum_b = calloc((size_t)b->n, sizeof(int));
    if (NULL == renum_b)
    {
      free(renum_a);
      assert(0);
      return FALSE;
    }
  }
  else
  {
    renum_a = NULL;
    renum_b = NULL;
  }

  status = rule_aggregate_cgraph_append_raw(a, rule_aggregate_full,
                                            b, rule_aggregate_b, c,
                                            renum_a, renum_b);

  if (NULL != renum_a)
  {
    free(renum_a);
  }

  if (NULL != renum_b)
  {
    free(renum_b);
  }

  return status;
}

/*
 * a                   : current full cgraph.
 * rule_aggregate_full : Current full rule aggregate.
 * b                   : Cgraph to merge.
 * rule_aggregate_b    : Rule aggregate associated with b.
 * c                   : Will contain fully merged cgraph.
 * renum_a             : Will map a nodes to c nodes.
 * renum_b             : Will map b nodes to c nodes.
 */
boolean
rule_aggregate_cgraph_append_raw(cgraph_t *a, int *rule_aggregate_full,
                                 cgraph_t *b, int *rule_aggregate_b,
                                 cgraph_t *c, int *renum_a, int *renum_b)
{
  int status;
  int *base_full;
  int *base_b;

  if (NULL != a)
  {
    /// Only append cgraphs if they exist.
    status = cgraph_append(a, renum_a, b, renum_b, c);
    if (FALSE == status)
    {
      return FALSE;
    }
  }

  /* Merge rule_aggregate_b into rule_aggregate_full, which will go with c. */
  base_full = rule_aggregate_full;
  base_b = rule_aggregate_b;

  // factor in multiplier.
  base_full[AGG_MULTIPLIER_OFFSET] *= base_b[AGG_MULTIPLIER_OFFSET];

  // Skip multiplier.
  base_full++;
  base_b++;

  // Split directives.
  status = rule_aggregate_double_cgraph_append(&base_full, &base_b,
                                               renum_a, renum_b);
  if (FALSE == status)
  {
    return FALSE;
  }

  // Bind directives.
  status = rule_aggregate_single_cgraph_append(&base_full, &base_b,
                                               renum_a, renum_b);
  if (FALSE == status)
  {
    return FALSE;
  }

  // Self bind directives.
  status = rule_aggregate_double_cgraph_append(&base_full, &base_b,
                                               renum_a, renum_b);
  if (FALSE == status)
  {
    return FALSE;
  }

  // State change directives.
  status = rule_aggregate_state_change_cgraph_append(&base_full, &base_b,
                                                     renum_a, renum_b);
  if (FALSE == status)
  {
    return FALSE;
  }

  // Transient bind directives.
  status = rule_aggregate_single_cgraph_append(&base_full, &base_b,
                                               renum_a, renum_b);
  if (FALSE == status)
  {
    return FALSE;
  }

  // Transient self bind directives.
  status = rule_aggregate_double_cgraph_append(&base_full, &base_b,
                                               renum_a, renum_b);
  if (FALSE == status)
  {
    return FALSE;
  }

  return TRUE;
}

static int
rule_aggregate_split(int *base_a, void *data1)
{
  cgraph_t *cg;
  int pri;
  int sec;

  cg = data1;

#ifndef DS_NO_DOUBLECHECK
  if((0 == base_a[AGG_PRI_USED_OFFSET]) ||
     (0 == base_a[AGG_SEC_USED_OFFSET]))
  {
    return 0;
  }
#endif

  pri = base_a[AGG_PRI_NODE_OFFSET];
  sec = base_a[AGG_SEC_NODE_OFFSET];

#ifndef DS_NO_DOUBLECHECK
  assert(pri != sec);
  assert(pri != -1);
  assert(sec != -1);
#endif

  // Break link between pri to sec and sec to pri.
  DELELEMENT(GRAPHROW(cg->g.a, pri, cg->m), sec);
  DELELEMENT(GRAPHROW(cg->g.a, sec, cg->m), pri);

  return 1;
}

static int
rule_aggregate_bind(int *base_a, void *data1)
{
  cgraph_t *cg;
  int pri;
  int sec;

  cg = data1;

#ifndef DS_NO_DOUBLECHECK
  assert((1 == base_a[AGG_PRI_USED_OFFSET]) &&
         (1 == base_a[AGG_SEC_USED_OFFSET]));
#endif

  pri = base_a[AGG_PRI_NODE_OFFSET];
  sec = base_a[AGG_SEC_NODE_OFFSET];

#ifndef DS_NO_DOUBLECHECK
  assert(pri != sec);
#endif

  // Create link between pri to sec and sec to pri.
#ifdef AGGRESSIVE_CHECK
assert(sec < cg->n);  
assert(pri < cg->n);  
#endif
  ADDELEMENT(GRAPHROW(cg->g.a, pri, cg->m), sec);
  ADDELEMENT(GRAPHROW(cg->g.a, sec, cg->m), pri);

  return 1;
}

static int
rule_aggregate_state_change(int *base_a, void *data1)
{
  cgraph_t *cg;
  int node;
  int old;
  int new;
  dynints_t *state_dynints;
  int state_index;

  cg = data1;

#ifndef DS_NO_DOUBLECHECK
  assert(1 == base_a[AGG_PRI_USED_OFFSET]);
#endif

  node = base_a[AGG_PRI_NODE_OFFSET];
  old = base_a[AGG_OLD_STATE_OFFSET];
  new = base_a[AGG_NEW_STATE_OFFSET];
  state_dynints = cg->state_dynints_array.a[node];

#ifndef DS_NO_DOUBLECHECK
  assert(old != new);
#endif

  // Find state that needs to be changed.
  state_index = dynints_int_find(old, state_dynints);

#ifndef DS_NO_DOUBLECHECK
  assert(state_index != -1);
#endif

  state_dynints->ints.a[state_index] = new;

  return 1;
}

boolean
rule_aggregate_cgraph_apply(cgraph_t *a, int *rule_aggregate_full)
{
  int *base;
  int status;

  base = rule_aggregate_full;

  // Skip multiplier.
  base++;

  // Split directives.
  status = rule_aggregate_iterate
    (rule_aggregate_split, &base, a, FULL_AGG_SPLIT_LENGTH);
  if (0 == status)
  {
    return 0;
  }

  // Bind directives.
  status = rule_aggregate_iterate
    (rule_aggregate_bind, &base, a, FULL_AGG_BIND_LENGTH);
  if (0 == status)
  {
    return 0;
  }

  // Self bind directives.
  status = rule_aggregate_iterate
    (rule_aggregate_bind, &base, a, FULL_AGG_BIND_LENGTH);
  if (0 == status)
  {
    return 0;
  }

  // State change directives.
  status = rule_aggregate_iterate
    (rule_aggregate_state_change, &base, a, FULL_AGG_STATE_CHANGE_LENGTH);
  if (0 == status)
  {
    return 0;
  }

  // Transient bind directives do not modify final graph.

  // Transient self bind directives do not modify final graph.


  return TRUE;
}

static boolean
rule_aggregate_raw_aggregate_compare
  (const int **pbase1, const int **pbase2,
   const size_t used_offset, const size_t node_offset,
   const size_t block_length)
{
  const int *base1;
  const int *base2;

  base1 = *pbase1;
  base2 = *pbase2;

  while(*base1 != -1)
  {
#ifndef DS_NO_DOUBLECHECK
    assert(*base2 != -1);
#endif

    if (base1[node_offset] != base2[node_offset])
    {
      return FALSE;
    }

    base1 += block_length;
    base2 += block_length;
  }

#ifndef DS_NO_DOUBLECHECK
  assert(-1 == *base1);
  assert(-1 == *base2);
#endif

  // Skip terminator.
  base1++;
  base2++;

  *pbase1 = base1;
  *pbase2 = base2;

  return TRUE;
}

static boolean
rule_aggregate_double_aggregate_compare
  (const int **pbase1, const int **pbase2)
{
  const int *base1;
  const int *base2;

  base1 = *pbase1;
  base2 = *pbase2;

  while(*base1 != -1)
  {
    int pri1 = base1[AGG_PRI_NODE_OFFSET];
    int pri2 = base2[AGG_PRI_NODE_OFFSET];
    int sec1 = base1[AGG_SEC_NODE_OFFSET];
    int sec2 = base2[AGG_SEC_NODE_OFFSET];

#ifndef DS_NO_DOUBLECHECK
    assert(*base2 != -1);
#endif

    if (!(((pri1 == pri2) && (sec1 == sec2)) ||
          ((pri1 == sec2) && (sec1 == pri2))))
    {
      return FALSE;
    }

    base1 += AGG_DOUBLE_LENGTH;
    base2 += AGG_DOUBLE_LENGTH;
  }

#ifndef DS_NO_DOUBLECHECK
  assert(-1 == *base1);
  assert(-1 == *base2);
#endif

  // Skip terminator.
  base1++;
  base2++;

  *pbase1 = base1;
  *pbase2 = base2;

  return TRUE;
}

static boolean
rule_aggregate_state_change_aggregate_compare
  (const int **pbase1, const int **pbase2)
{
  int status;

  status = rule_aggregate_raw_aggregate_compare
             (pbase1, pbase2,
              AGG_PRI_USED_OFFSET, AGG_PRI_NODE_OFFSET,
              AGG_STATE_CHANGE_LENGTH);

  if (FALSE == status)
  {
    return FALSE;
  }

  return TRUE;
}

static boolean
rule_aggregate_single_aggregate_compare
  (const int **pbase1, const int **pbase2)
{
  int status;

  status = rule_aggregate_raw_aggregate_compare
             (pbase1, pbase2,
              AGG_PRI_USED_OFFSET, AGG_PRI_NODE_OFFSET,
              AGG_SINGLE_LENGTH);

  if (FALSE == status)
  {
    return FALSE;
  }

  return TRUE;
}

/*
** Returns FALSE if any part of the aggregate is different.
** For split or self bind, ab is treated same as ba to properly
** handle powers of two.
*/
boolean
rule_aggregate_compare(const int *aggregate1, const int *aggregate2)
{
  int status;
  const int *base1;
  const int *base2;

  // Skip multiplier.
  base1 = aggregate1 + 1;
  base2 = aggregate2 + 1;

  // Split directives.
  status = rule_aggregate_double_aggregate_compare
             (&base1, &base2);
  if (FALSE == status)
  {
    return FALSE;
  }

  // Bind directives.
  status = rule_aggregate_single_aggregate_compare
             (&base1, &base2);
  if (FALSE == status)
  {
    return FALSE;
  }

  // Self bind directives.
  status = rule_aggregate_double_aggregate_compare
             (&base1, &base2);
  if (FALSE == status)
  {
    return FALSE;
  }

  // State change directives.
  status = rule_aggregate_state_change_aggregate_compare
             (&base1, &base2);
  if (FALSE == status)
  {
    return FALSE;
  }

  // Transient bind directives.
  status = rule_aggregate_single_aggregate_compare
             (&base1, &base2);
  if (FALSE == status)
  {
    return FALSE;
  }

  // Transient self bind directives.
  status = rule_aggregate_double_aggregate_compare
             (&base1, &base2);
  if (FALSE == status)
  {
    return FALSE;
  }

  return TRUE;
}

static int *
rule_aggregate_output_single(const output_type ptype, const int level,
                             int *base)
{
  // Print start.
  CUSTOM_PRINTF(ptype, level, "( ");

  while (*base != -1)
  {
    if (0 != base[AGG_USED_OFFSET])
    {
      CUSTOM_PRINTF(ptype, level, "%d ",
                    base[AGG_USED_OFFSET] ? base[AGG_NODE_OFFSET] : -1);
    }
    else
    {
      CUSTOM_PRINTF(ptype, level, "? ");
    }
    base += AGG_SINGLE_LENGTH;
  }

  // Print end.
  CUSTOM_PRINTF(ptype, level, ")");

  // Skip terminator.
  base++;

  return base;
}

static int *
rule_aggregate_output_double(const output_type ptype, const int level,
                             int *base)
{
  // Print start.
  CUSTOM_PRINTF(ptype, level, "( ");

  while (*base != -1)
  {
    if (0 != base[AGG_PRI_USED_OFFSET])
    {
      CUSTOM_PRINTF(ptype, level, "%d:", base[AGG_PRI_NODE_OFFSET]);
    }
    else
    {
      CUSTOM_PRINTF(ptype, level, "?:");
    }

    if (0 != base[AGG_SEC_USED_OFFSET])
    {
      CUSTOM_PRINTF(ptype, level, "%d ", base[AGG_SEC_NODE_OFFSET]);
    }
    else
    {
      CUSTOM_PRINTF(ptype, level, "? ");
    }

    base += AGG_DOUBLE_LENGTH;
  }

  // Print end.
  CUSTOM_PRINTF(ptype, level, ")");

  // Skip terminator.
  base++;

  return base;
}

static int *
rule_aggregate_output_state_change(const output_type ptype, const int level,
                                   int *base)
{
  // Print start
  CUSTOM_PRINTF(ptype, level, "( ");

  while (*base != -1)
  {
    if (TRUE == base[AGG_USED_OFFSET])
    {
      CUSTOM_PRINTF(ptype, level, "%d:%d->%d ", base[AGG_NODE_OFFSET],
                    base[AGG_OLD_STATE_OFFSET], base[AGG_NEW_STATE_OFFSET]);
    }
    else
    {
      CUSTOM_PRINTF(ptype, level, "?:?->? ");
    }
    base += AGG_STATE_CHANGE_LENGTH;
  }

  // Print end.
  CUSTOM_PRINTF(ptype, level, ")");

  // Skip terminator.
  base++;

  return base;
}

void
rule_aggregate_output_raw(const output_type ptype, const int level,
                          int *rule_aggregate, size_t num_ints)
{
  size_t i;
  int *pcurrent_int;

  for (i = 0, pcurrent_int = rule_aggregate;
       i < num_ints;
       i++, pcurrent_int++)
  {
    CUSTOM_PRINTF(ptype, level, "%d ", *pcurrent_int);
  }
}

void
rule_aggregate_output(const output_type ptype, const int level,
                      int *rule_aggregate)
{
  int *base;

  base = rule_aggregate;

  // Output multiplier.
  CUSTOM_PRINTF(ptype, level, "%d:", base[AGG_MULTIPLIER_OFFSET]);
  base++;

  // Output split block.
  base = rule_aggregate_output_double(ptype, level, base);

  // Output bind block.
  base = rule_aggregate_output_single(ptype, level, base);

  // Output self bind block.
  base = rule_aggregate_output_double(ptype, level, base);

  // Output state change block.
  base = rule_aggregate_output_state_change(ptype, level, base);

  // Output transient bind block.
  base = rule_aggregate_output_single(ptype, level, base);

  // Output transient self bind block.
  base = rule_aggregate_output_double(ptype, level, base);
}

void
rule_aggregate_output_full(const output_type ptype, const int level,
                           int *rule_aggregate)
{
  int *base;

  base = rule_aggregate;

  // Output multiplier.
  CUSTOM_PRINTF(ptype, level, "%d:", base[AGG_MULTIPLIER_OFFSET]);
  base++;

  // Output split block.
  base = rule_aggregate_output_double(ptype, level, base);

  // Output bind block.
  base = rule_aggregate_output_double(ptype, level, base);

  // Output self bind block.
  base = rule_aggregate_output_double(ptype, level, base);

  // Output state change block.
  base = rule_aggregate_output_state_change(ptype, level, base);

  // Output transient bind block.
  base = rule_aggregate_output_double(ptype, level, base);

  // Output transient self bind block.
  base = rule_aggregate_output_double(ptype, level, base);
}

int
rule_aggregate_calculate_product_count(cgraph_t *a, int *rule_aggregate)
{
  int product_count;
  int *base;
  int status;

  // Perform operations on specie graph.
  base = rule_aggregate;

  // Skip multiplier.
  base++;

  // Split directives.
  status = rule_aggregate_iterate
    (rule_aggregate_split, &base, a, AGG_SPLIT_LENGTH);
  if (0 == status)
  {
    return 0;
  }

  // Skip bind directives.
  rule_aggregate_skip(&base, AGG_BIND_LENGTH);

  // Self bind directives.
  status = rule_aggregate_iterate
    (rule_aggregate_bind, &base, a, AGG_SELF_BIND_LENGTH);
  if (0 == status)
  {
    return 0;
  }


  // Check if graph is split or not.
  product_count = cgraph_seperate(a, NULL, NULL);


  // Restore specie graph.
  base = rule_aggregate;

  // Skip multiplier.
  base++;

  // Reverse split directives.
  status = rule_aggregate_iterate
    (rule_aggregate_bind, &base, a, AGG_SPLIT_LENGTH);
  if (0 == status)
  {
    return 0;
  }

  // Skip bind directives.
  rule_aggregate_skip(&base, AGG_BIND_LENGTH);

  // Reverse bind directives.
  status = rule_aggregate_iterate
    (rule_aggregate_split, &base, a, AGG_SELF_BIND_LENGTH);
  if (0 == status)
  {
    return 0;
  }

  return product_count;
}

