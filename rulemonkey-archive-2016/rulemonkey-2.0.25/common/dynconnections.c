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
#include "dynconnections.h"
#include "output.h"

void suppress_warnings28() {SUPPRESS_WARNINGS;}

// Local function declarations.
static boolean dynconnections_resize(dynconnections_t *dynconnections,
                                     const size_t count);
static boolean dynconnections_connection_exists
         (dynconnections_t *dynconnections, const int a, const int b);
static int dynconnections_same_aggregate(const size_t rule_index_a,
                                         const size_t rule_index_b,
                                         int index_count,
                                         const int *specie_rule_indexes);

/*
** Used when building specie graphs and completing rules.
*/
boolean
dynconnections_connection_add(const int specie_rule_index, const char *name,
                              const int node, dynconnections_t *dynconnections)
{
  size_t i;
  boolean status;

  // Get index of existing, or create new entry.
  i = dynconnections_get_index(specie_rule_index, name, dynconnections);
  if (0 == i)
  {
    assert(0);
    return FALSE;
  }

  status = dynints_int_add(node, dynconnections->dynints.a[i]);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

  if (dynconnections->dynints.a[i]->count > 2)
  {
    error_printf("Found more than two instances of connection \"%s\"\n", name);
    return FALSE;
  }

  return TRUE;
}

void
dynconnections_output(const output_type ptype, const int level,
                      dynconnections_t *dynconnections)
{
  dynints_t **dynints;
  int i;
  int j;

  dynints = dynconnections->dynints.a;

  CUSTOM_PRINTF(ptype, level, "Dynconnections\n");
  for (i = 1; i < dynconnections->count; i++)
  {
    size_t ints_count;
    int *ints;

    ints_count = dynints[i]->count;
    ints = dynints[i]->ints.a;

    CUSTOM_PRINTF(ptype, level, "  %d-%s:",
                  dynconnections->specie_rule_indexes.a[i],
                  dynconnections->names.a[i]);
    for (j = 0; j < ints_count; j++)
    {
      CUSTOM_PRINTF(ptype, level, " %d", ints[j]);
    }
    CUSTOM_PRINTF(ptype, level, "\n");
  }
}

boolean
dynconnections_compare_rules(size_t lhs_rule_index, size_t rhs_rule_index,
                             size_t rhs_index,
                             dynconnections_t *dynconnections)
{
  size_t i;
  size_t count;
  dynints_t **pcurrent_dynint;
  int lhs_single;
  int rhs_single;

  lhs_single = 0;
  rhs_single = 0;

  count = dynconnections->count;
  // Index 0 is reserved for errors.
  for (i = 1, pcurrent_dynint = dynconnections->dynints.a + 1;
       i < count;
       i++, pcurrent_dynint++)
  {
    size_t ints_count;
    int *ints;

    ints_count = (*pcurrent_dynint)->count;
    ints = (*pcurrent_dynint)->ints.a;

    if (1 == ints_count)
    {
      if (ints[0] == lhs_rule_index)
      {
        // Solitary connection is on left hand side.
        lhs_single = 1;
      }
      else if (ints[0] == rhs_rule_index)
      {
        // Solitary connection is on left hand side.
        rhs_single = 1;
      }
    }
  }

  if (rhs_single != lhs_single)
  {
    // Single connections must match.
    return FALSE;
  }

  // No unshared connection found.
  return TRUE;
}

void
dynconnections_apply_to_graph(cgraph_t *cg, dynconnections_t *dynconnections)
{
  size_t i;
  size_t count;
  dynints_t **dynints;

  count = dynconnections->count;
  dynints = dynconnections->dynints.a;
  // 0 is reserved to signify an error.
  for (i = 1; i < count; i++)
  {
    size_t ints_count;
    size_t j;
    int *ints;

    ints_count = dynints[i]->count;
    ints = dynints[i]->ints.a;

    /* Fully connect all nodes listed. */
    for (j = 0; j < ints_count; j++)
    {
      size_t k;

      for (k = 0; k < ints_count; k++)
      {
        /* Don't connect node to itself. */
        if (j != k)
        {
          cgraph_connect_nodes(cg, ints[j], ints[k]);
        }
      }
    }
  }
}

/**
 * Apply connection-connection bindings.
 */
boolean
dynconnections_apply_to_rules(const int lhs_count,
                              const int rhs_count,
                              int strict,
                              const int *specie_rule_indexes,
                              dyncolors_t *dyncolors,
                              dynrules_t *dynrules,
                              dynconnections_t *dynconnections)
{
  size_t i;
  size_t j;
  size_t count;
  size_t first_rule_index;
  size_t last_rule_index;
  size_t rhs_rule_index;
  dynints_t **dynints;
  

  rhs_rule_index = specie_rule_indexes[lhs_count];
  first_rule_index = specie_rule_indexes[0];
  last_rule_index = dynrules->count - 1;

  count = dynconnections->count;
  dynints = dynconnections->dynints.a;

  // 0 is reserved to signify an error.
  // Apply all connections.
  for (i = 1; i < count; i++)
  {
    size_t ints_count;
    int *ints;
    const char *name;

    ints_count = dynints[i]->count;
    ints = dynints[i]->ints.a;
    name = dynconnections->names.a[i];

    // Sanity check.
    for (j = 0; j < ints_count; j++)
    {
      if ((ints[j] < first_rule_index) || (ints[j] > last_rule_index))
      {
        // Rule out of range.
        assert(0);
        return FALSE;
      }
    }

    if ((0 == strcmp(name, "+")) ||
        ((FALSE == strict) && (1 == ints_count)))
    {
      int m;

      // Connections that can be applied all at once on left and right.
      // Specify one or more connections present
      for (m = 0; m < ints_count; m++)
      {
        boolean status;

        status = dynrules_rule_connect_any_one_rule
                   (ints[m], dynrules);
        if (FALSE == status)
        {
          assert(0);
          return FALSE;
        }
      }
    }
    else if (0 == strcmp(name, "?"))
    {
      int m;

      // Connections that can be applied all at once on left and right.
      // Specify zero or one connection present
      for (m = 0; m < ints_count; m++)
      {
        boolean status;

        status = dynrules_rule_connect_zero_or_one_rule
                   (ints[m], dynrules);
        if (FALSE == status)
        {
          assert(0);
          return FALSE;
        }
      }
    }
    else if ((2 == ints_count))
    {
      // Apply pair.
      if ((ints[0] < rhs_rule_index) && (ints[1] < rhs_rule_index))
      {
        // Pair on left side.
        // Link rules together.
        dynrules_rule_add_adjacent_rule(ints[0], ints[1], dynrules);
        dynrules_rule_add_adjacent_rule(ints[1], ints[0], dynrules);
      }
      else if ((ints[0] >= rhs_rule_index) && (ints[1] >= rhs_rule_index))
      {
        // Pair on right side.
        // Link rules together.
        dynrules_rule_add_adjacent_rule(ints[0], ints[1], dynrules);
        dynrules_rule_add_adjacent_rule(ints[1], ints[0], dynrules);
      }
      else
      {
        // Pair not on the same side of reaction.
        return FALSE;
      }
    }
    else
    {
      // Single connection specifier specifies wildcard connection.
    }
  }

  return TRUE;
}

boolean
dynconnections_apply_map(const int lhs_count, const int rhs_count,
                         const int *specie_rule_indexes,
                         dynrules_t *dynrules, dynmaps_t *dynmaps,
                         dynconnections_t *dynconnections)
{
  size_t i;
  size_t count;
  const char *name;
  dynints_t **dynints;
  int *psplit_index;
  int f_split_index;
  int r_split_index;
  int *pbind_index;
  int f_bind_index;
  int r_bind_index;
  int *pself_bind_index;
  int f_self_bind_index;
  int r_self_bind_index;
  size_t rhs_index;
#if 0
  /* Need to finish transient_* support. */
  int *transient_bind_index
  int f_transient_bind_index;
  int r_transient_bind_index;
  int *transient_self_bind_index;
  int f_transient_self_bind_index;
  int r_transient_self_bind_index;
#endif

  f_split_index = 0;
  r_split_index = 0;
  f_bind_index = 0;
  r_bind_index = 0;
  f_self_bind_index = 0;
  r_self_bind_index = 0;
#if 0
  f_transient_bind_index = 0;
  r_transient_bind_index = 0;
  f_transient_self_bind_index = 0;
  r_transient_self_bind_index = 0;
#endif

  rhs_index = (unsigned)specie_rule_indexes[lhs_count];

  count = dynconnections->count;
  dynints = dynconnections->dynints.a;
  for (i = 1, dynints = dynconnections->dynints.a + 1;
       i < count;
       i++, dynints++)
  {
    int status;
    size_t ints_count;
    int *ints;
    size_t rule_index_a;
    size_t rule_index_b;

    ints_count = (*dynints)->count;
    ints = (*dynints)->ints.a;
    name = dynconnections->names.a[i];

    if ((0 == strcmp(name, "+")) || (0 == strcmp(name, "?")))
    {
      // Nothing to do here, move along.
      continue;
    }

    // Check if link is present on other side of reaction.

    if (2 == ints_count)
    {
      if ((ints[0] < rhs_index) && (ints[1] < rhs_index))
      {
        status = dynmaps_dest_from_source(&rule_index_a, (size_t)ints[0],
                                          dynmaps);
        if (FALSE == status)
        {
          assert(0);
          return FALSE;
        }
        status = dynmaps_dest_from_source(&rule_index_b, (size_t)ints[1],
                                          dynmaps);
        if (FALSE == status)
        {
          assert(0);
          return FALSE;
        }

        // Check if components on other side of reaction are already connected.
        if (1 == dynconnections_connection_exists(dynconnections,
                                                  rule_index_a, rule_index_b))
        {
          // Link exists on both sides of reaction, so nothing to do.
          continue;
        }
        else
        {
          // Link only on left side of reaction.  Need to add split directive
          // to molecules with link that will be removed, and bind directive
          // to matching molecules that are missing the link.

          // Existing link on lhs, so will add bind directives to rhs.
          psplit_index = &f_split_index;
          pbind_index = &r_bind_index;
          pself_bind_index = &r_self_bind_index;
        }
      }
      else if ((ints[0] >= rhs_index) && (ints[1] >= rhs_index))
      {
        status = dynmaps_source_from_dest(&rule_index_a, (size_t)ints[0],
                                          dynmaps);
        if (FALSE == status)
        {
          return FALSE;
        }
        status = dynmaps_source_from_dest(&rule_index_b, (size_t)ints[1],
                                          dynmaps);
        if (FALSE == status)
        {
          return FALSE;
        }

        if (1 == dynconnections_connection_exists(dynconnections,
                                                  rule_index_a, rule_index_b))
        {
          // Link exists on both sides of reaction, so nothing to do.
          continue;
        }
        else
        {
          // Link only on right side of reaction.  Need to add split directive
          // to molecules with link that will be removed, and bind directive
          // to matching molecules that are missing the link.

          // Existing link on rhs, so will add bind directives to rhs.
          psplit_index = &r_split_index;
          pbind_index = &f_bind_index;
          pself_bind_index = &f_self_bind_index;
        }
      }
      else
      {
        // Single bind point on each side of reaction simply used to designate
        // something present, not to specify a full connection.
        continue;
      }

      // At this point, ints[0] & ints[1] are bound, but will
      // be split when rule is applied.
      // Also, rule_index_a & rule_index_b are not bound yet, but will be
      // bound when rule is applied.

      // Note split points in rules.
      status = dynrules_rule_add_primary_split_index(ints[0], *psplit_index,
                                                     dynrules);
      if (FALSE == status)
      {
        /* Error specifying split point in rules. */
        assert(0);
        return FALSE;
      }
      status = dynrules_rule_add_secondary_split_index(ints[1], *psplit_index,
                                                       dynrules);
      if (FALSE == status)
      {
        /* Error specifying split point in rules. */
        assert(0);
        return FALSE;
      }

      *psplit_index += 1;

      if (1 == dynconnections_same_aggregate(rule_index_a, rule_index_b,
                                             lhs_count + rhs_count,
                                             specie_rule_indexes))
      {
        // Both points are on same entity, so it will form a ring.

        // Note bind points in rules.
        status = dynrules_rule_add_primary_self_bind_index
                   ((int)rule_index_a, *pself_bind_index, dynrules);
        if (FALSE == status)
        {
          assert(0);
          return FALSE;
        }
        status = dynrules_rule_add_secondary_self_bind_index
                   ((int)rule_index_b, *pself_bind_index, dynrules);
        if (FALSE == status)
        {
          assert(0);
          return FALSE;
        }

        *pself_bind_index += 1;
      }
      else
      {
        // Bind points are on different entities, so will form aggregate.

        // Note bind points in rules.
        status = dynrules_rule_add_bind_index((int)rule_index_a,
                                              *pbind_index, dynrules);
        if (FALSE == status)
        {
          assert(0);
          return FALSE;
        }
        status = dynrules_rule_add_bind_index((int)rule_index_b,
                                              *pbind_index, dynrules);
        if (FALSE == status)
        {
          assert(0);
          return FALSE;
        }

        *pbind_index += 1;
      }
    }
  }

  return TRUE;
}

/*
** Check if rule_a and rule_b are part of the same aggregate or not.
*/
static int
dynconnections_same_aggregate(const size_t rule_index_a,
                              const size_t rule_index_b,
                              int index_count,
                              const int *specie_rule_indexes)
{
  size_t index_x;
  size_t index_y;
  int current_index;


  // Make sure index_x is less than index_y.
  if (rule_index_a < rule_index_b)
  {
    index_x = rule_index_a;
    index_y = rule_index_b;
  }
  else
  {
    index_x = rule_index_b;
    index_y = rule_index_a;
  }

  current_index = 0;
  while ((current_index < index_count) &&
         (specie_rule_indexes[current_index] < index_x))
  {
    current_index++;
  }

  if (current_index == index_count)
  {
    // Both rules are part of last aggregate.
    return 1;
  }

  /*
  ** At this point, index_x is part of aggregate starting at
  ** specie_rule_indexes[current_index - 1] and ending at
  ** specie_rule_indexes[current_index].
  */
  if (specie_rule_indexes[current_index] > index_y)
  {
    // index_y is in same aggregate as index_x.
    return 1;
  }

  return 0;
}

size_t
dynconnections_connection_alloc(dynconnections_t *dynconnections)
{
  int status;

  status = dynconnections_resize(dynconnections, dynconnections->count + 1);
  if (FALSE == status)
  {
    assert(0);
    return 0;
  }

  dynconnections->dynints.a[dynconnections->count] =
    dynints_create(INITIAL_DYNINTS);
  if (NULL == dynconnections->dynints.a[dynconnections->count])
  {
    assert(0);
    return 0;
  }

  return dynconnections->count++;
}

size_t
dynconnections_get_index(const int specie_rule_index, const char *name,
                         dynconnections_t *dynconnections)
{
  size_t i;
  size_t count;
  char **pcurrent_name; 
  int *pcurrent_specie_rule_index;

  if (NULL == name)
  {
    assert(0 && "dynconnections_get_index got an empty name pointer");
    return 0;
  }

  count = dynconnections->count;
  // Index 0 is reserved for errors.
  for (i = 1, pcurrent_name = dynconnections->names.a + 1,
       pcurrent_specie_rule_index = dynconnections->specie_rule_indexes.a + 1;
       i < count;
       i++, pcurrent_name++, pcurrent_specie_rule_index++)
  {   
    if ((specie_rule_index == *pcurrent_specie_rule_index) &&
        (0 == strcmp(name, *pcurrent_name)))
    {
      // Found a match.
      return i;
    }
  }
  
  // Name was not found, need to create new connection.
  i = dynconnections_connection_alloc(dynconnections);
  if (0 == i)
  {
    return 0;
  }

  dynconnections->specie_rule_indexes.a[i] = specie_rule_index;
  dynconnections->names.a[i] = strdup(name);
  if (NULL == dynconnections->names.a[i])
  {
    return 0;
  }

  return i;
}

dynconnections_t *
dynconnections_create(const size_t input_init_count)
{
  dynconnections_t *new_dynconnections;
  int status;
  size_t init_count;

  if (0 == input_init_count)
  {
    // Index 0 is reserved for errors, make sure it exists.
    init_count = 1;
  }
  else
  {
    init_count = input_init_count;
  }

  // All fields need to be initialized to zero anyways, so use calloc.
  new_dynconnections = calloc(1, sizeof(dynconnections_t));
  if (NULL == new_dynconnections)
  {
    assert(0);
    error_printf("Unable to allocate new dynconnections.\n");
    return NULL;
  }

  status =
    dynarray_create3(init_count,
                     &new_dynconnections->alloc_count,
                     sizeof(int),
                     &new_dynconnections->specie_rule_indexes.v,
                     sizeof(char *), &new_dynconnections->names.v,
                     sizeof(dynints_t *),
                     &new_dynconnections->dynints.v);
  if (0 == status)
  {
    // Error message output in dynarray_create2.
    assert(0);
    free(new_dynconnections);
    return NULL;
  }

  // Reserve 0 index for errors.
  new_dynconnections->count = 1;
  new_dynconnections->specie_rule_indexes.a[0] = -1;
  new_dynconnections->names.a[0] = "myerror";
  new_dynconnections->dynints.a[0] = NULL;

  return new_dynconnections;
}

static boolean
dynconnections_connection_exists(dynconnections_t *dynconnections,
                                 const int a, const int b)
{
  size_t i;
  size_t count;
  char **pcurrent_name; 
  dynints_t **pdynints;

  count = dynconnections->count;
  // Index 0 is reserved for errors.
  for (i = 1, pcurrent_name = dynconnections->names.a + 1,
       pdynints = dynconnections->dynints.a + 1;
       i < count;
       i++, pcurrent_name++, pdynints++)
  {   
    dynints_t *dynints;
    int *ints;


    dynints = *pdynints;
    ints = dynints->ints.a;

    if ((2 == dynints->count) &&
        (((ints[0] == a) && (ints[1] == b)) ||
         ((ints[0] == b) && (ints[1] == a))))
    {
      // Found a match.
      return TRUE;
    }
  }
  
  // Match not found.
  return FALSE;
}

static boolean
dynconnections_resize(dynconnections_t *dynconnections, const size_t count)
{
  int status;

  status =
    dynarray_resize3(count,
                     &dynconnections->alloc_count,
                     sizeof(int),
                     &dynconnections->specie_rule_indexes.v,
                     sizeof(char *), &dynconnections->names.v,
                     sizeof(dynints_t *), &dynconnections->dynints.v);
  if (0 == status)
  {
    // Error message output in dynarray_resize2.
    return 0;
  }

  return 1;
}

void
dynconnections_reset(dynconnections_t *dynconnections)
{
  size_t i;
  int *specie_rule_indexes;
  char **names;
  dynints_t **dynints;
  size_t count;

  specie_rule_indexes = dynconnections->specie_rule_indexes.a;
  names = dynconnections->names.a;
  dynints = dynconnections->dynints.a;
  count = dynconnections->count;
  // Index 0 is reserved for errors.
  for (i = 1; i < count; i++)
  {
    free(names[i]);
    specie_rule_indexes[i] = -1;
    names[i] = NULL;
    dynints_destroy(dynints[i]);
  }

  dynconnections->count = 1;
}

void
dynconnections_destroy(dynconnections_t *dynconnections)
{
  size_t i;
  char **names;
  dynints_t **dynints;
  size_t count;

  names = dynconnections->names.a;
  dynints = dynconnections->dynints.a;
  count = dynconnections->count;
  // Index 0 is reserved for errors.
  for (i = 1; i < count; i++)
  {
    free(names[i]);
    dynints_destroy(dynints[i]);
  }

  free(dynconnections->specie_rule_indexes.a);
  free(dynconnections->names.a);
  free(dynconnections->dynints.a);
  free(dynconnections);
}

