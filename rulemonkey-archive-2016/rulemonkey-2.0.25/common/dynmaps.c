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
#include "dynmaps.h"
#include "dynarray.h"
#include "constants.h"
#include "output.h"
#include "dynrules.h"

void suppress_warnings30() {SUPPRESS_WARNINGS;}

// Local function declarations.
static boolean dynmaps_resize(dynmaps_t *dynmaps, const size_t count);

boolean
dynmaps_map_add(const size_t source, const size_t dest, dynmaps_t *dynmaps)
{
  boolean status;
  size_t new_index;

  status = dynmaps_map_alloc(&new_index, dynmaps);
  if (FALSE == status)
  {
    return FALSE;
  }

  dynmaps->sources.a[new_index] = source;
  dynmaps->destinations.a[new_index] = dest;

  return TRUE;
}

void
dynmaps_output(const output_type ptype, const int level, dynmaps_t *dynmaps)
{
  size_t i;

  CUSTOM_PRINTF(ptype, level, "Dynmaps:\n");
  for (i = 0; i < dynmaps->count; i++)
  {
    CUSTOM_PRINTF(ptype, level, "  %lu <-> %lu\n",
                  dynmaps->sources.a[i], dynmaps->destinations.a[i]);
  }
}

/*
** dynmaps_apply_to_rules
**
** Add state changes to rules.
**
*/
boolean
dynmaps_apply_to_rules(dynmaps_t *dynmaps, dynrules_t *dynrules)
{
  size_t i;
  size_t j;
  size_t count;
  size_t *psource;
  size_t *pdest;
  int state_change_index;

  state_change_index = 0;

  count = dynmaps->count;

  for (i = 0, psource = dynmaps->sources.a, pdest = dynmaps->destinations.a;
       i < count;
       i++, psource++, pdest++)
  {
    rule_t *source_rule;
    rule_t *dest_rule;
    dynints_t *source_state_dynints;
    size_t source_state_count;
    int *psource_state_int;
    dynints_t *dest_state_dynints;
    size_t dest_state_count;
    int *pdest_state_int;

    source_rule = dynrules->rules.a + *psource;
    dest_rule = dynrules->rules.a + *pdest;

    source_state_dynints = source_rule->state_dynints;
    dest_state_dynints = dest_rule->state_dynints;

    source_state_count = source_state_dynints->count;
    dest_state_count = dest_state_dynints->count;
    if (source_state_count != dest_state_count)
    {
      // Mapped rules must always have the same number of states.
      return FALSE;
    }

    source_state_dynints = source_rule->state_dynints;
    dest_state_dynints = dest_rule->state_dynints;
    for (j = 0,
           psource_state_int = source_state_dynints->ints.a,
           pdest_state_int = dest_state_dynints->ints.a;
         j < source_state_count;
         j++, psource_state_int++, pdest_state_int++)
    {
      if (*psource_state_int != *pdest_state_int)
      {
        // State has changed.
        rule_add_modify_state(state_change_index, *psource_state_int,
                              *pdest_state_int, source_rule);
        rule_add_modify_state(state_change_index, *pdest_state_int,
                              *psource_state_int, dest_rule);
        state_change_index += 1;
      }
    }
  }

  return TRUE;
}

boolean
dynmaps_dest_from_source(size_t *pdest, const size_t source,
                         dynmaps_t *dynmaps)
{
  size_t i;
  size_t count;
  size_t *current_source;

  count = dynmaps->count;
  for (i = 0, current_source = dynmaps->sources.a;
       i < count;
       i++, current_source++)
  {
    if (*current_source == source)
    {
      // Found.
      *pdest = dynmaps->destinations.a[i];
      return TRUE;
    }
  }

  // Not found.
  assert(0);
  return FALSE;
}
boolean
dynmaps_source_from_dest(size_t *psource, const size_t dest,
                         dynmaps_t *dynmaps)
{
  size_t i;
  size_t count;
  size_t *current_dest;

  count = dynmaps->count;
  for (i = 0, current_dest = dynmaps->destinations.a;
       i < count;
       i++, current_dest++)
  {
    if (*current_dest == dest)
    {
      // Found.
      *psource = dynmaps->sources.a[i];
      return TRUE;
    }
  }

  // Not found.
  return FALSE;
} 

boolean
dynmaps_dest_exists(const size_t dest, dynmaps_t *dynmaps)
{
  size_t i;
  size_t count;
  size_t *current_dest;

  count = dynmaps->count;
  current_dest = dynmaps->destinations.a;
  for (i = 0, current_dest = dynmaps->destinations.a;
       i < count;
       i++, current_dest++)
  {
    if (*current_dest == dest)
    {
      // Found.
      return TRUE;
    }
  }

  // Not found.
  return FALSE;
}

boolean
dynmaps_map_alloc(size_t *pnew_index, dynmaps_t *dynmaps)
{
  int status;

  status = dynmaps_resize(dynmaps, dynmaps->count + 1);
  if (0 == status)
  {
    assert(0);
    return FALSE;
  }

  *pnew_index = dynmaps->count++;

  return TRUE;
}

dynmaps_t *
dynmaps_create(const size_t init_count)
{
  dynmaps_t *new_dynmaps;
  int status;

  new_dynmaps = calloc(1, sizeof(dynmaps_t));
  if (NULL == new_dynmaps)
  {
    error_printf("Unable to allocate new dynmaps.\n");
    return NULL;
  }

  status =
    dynarray_create2(init_count, &new_dynmaps->alloc_count,
                     sizeof(size_t), &new_dynmaps->sources.v,
                     sizeof(size_t), &new_dynmaps->destinations.v);
  if (0 == status)
  {
    // Error message output in dynarray_create.
    free(new_dynmaps);
    return NULL;
  }

  new_dynmaps->count = 0;

  return new_dynmaps;
}

static boolean
dynmaps_resize(dynmaps_t *dynmaps, const size_t count)
{
  int status;

  status = dynarray_resize2(count, &dynmaps->alloc_count,
                            sizeof(size_t), &dynmaps->sources.v,
                            sizeof(size_t), &dynmaps->destinations.v);
  if (0 == status)
  {
    // Error message output in dynarray_resize.
    return FALSE;
  }

  return TRUE;
}

void
dynmaps_destroy(dynmaps_t *dynmaps)
{
  free(dynmaps->sources.a);
  free(dynmaps->destinations.a);
  free(dynmaps);
}
