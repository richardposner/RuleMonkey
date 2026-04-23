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
#include "dynps_reactant_containers.h"
#include "dynarray.h"
#include "constants.h"
#include "rule_aggregate.h"
#include "output.h"

void suppress_warnings20() {SUPPRESS_WARNINGS;}

static boolean
dynps_reactant_containers_resize(dynps_reactant_containers_t *dynps_reactant_containers, const size_t count);

boolean
dynps_reactant_containers_ps_reactant_container_add
  (reaction_definition_t *reaction_definition, const int reactant_index,
   const int primary_particle, const size_t aggregate_size,
   const size_t full_aggregate_size, const int *aggregate,
   dynps_reactant_containers_t *dynps_reactant_containers)
{
  ps_reactant_container_t *new_ps_reactant_container;

  new_ps_reactant_container =
    dynps_reactant_containers_ps_reactant_container_alloc(dynps_reactant_containers);
  if (NULL == new_ps_reactant_container)
  {
    return FALSE;
  }

  ps_reactant_container_setup(new_ps_reactant_container, reaction_definition,
                              reactant_index, primary_particle, aggregate_size,
                              full_aggregate_size,
                              aggregate);

  return TRUE;
}

ps_reactant_container_t *
dynps_reactant_containers_ps_reactant_container_alloc
  (dynps_reactant_containers_t *dynps_reactant_containers)
{
  int status;

  status =
    dynps_reactant_containers_resize
      (dynps_reactant_containers, dynps_reactant_containers->count + 1);
  if (0 == status)
  {
    assert(0);
    return NULL;
  }

  return &dynps_reactant_containers->
           ps_reactant_containers.a[dynps_reactant_containers->count++];
}

dynps_reactant_containers_t *
dynps_reactant_containers_create(const size_t init_count)
{
  dynps_reactant_containers_t *new_dynps_reactant_containers;
  int status;

  new_dynps_reactant_containers =
    calloc(1, sizeof(dynps_reactant_containers_t));
  if (NULL == new_dynps_reactant_containers)
  {
    error_printf("Unable to allocate new dynps_reactant_containers.\n");
    return NULL;
  }

  status =
    dynarray_create
      (init_count, &new_dynps_reactant_containers->alloc_count,
       sizeof(ps_reactant_container_t),
       &new_dynps_reactant_containers->ps_reactant_containers.v);
  if (0 == status)
  {
    // Error message output in dynarray_create.
    free(new_dynps_reactant_containers);
    return NULL;
  }

  new_dynps_reactant_containers->count = 0;

  return new_dynps_reactant_containers;
}

/** @brief Reset object so that it can be reused later.
*/
void
dynps_reactant_containers_reset
  (dynps_reactant_containers_t *dynps_reactant_containers)
{
  size_t i;
  ps_reactant_container_t *current_rc;

  // Loop through alloc_count because containers may have been
  // initialized previously, but never reused.
  for (i = 0, current_rc = dynps_reactant_containers->ps_reactant_containers.a;
       i < dynps_reactant_containers->alloc_count;
       i++, current_rc++)
  {
    ps_reactant_container_destroy(current_rc);
  }

  dynps_reactant_containers->count = 0;
}

static boolean
dynps_reactant_containers_resize
  (dynps_reactant_containers_t *dynps_reactant_containers, const size_t count)
{
  int status;

  status =
    dynarray_resize
      (count, &dynps_reactant_containers->alloc_count,
       sizeof(ps_reactant_container_t),
       &dynps_reactant_containers->ps_reactant_containers.v);
  if (0 == status)
  {
    // Error message output in dynarray_resize.
    return FALSE;
  }

  return TRUE;
}

void
dynps_reactant_containers_destroy
  (dynps_reactant_containers_t *dynps_reactant_containers)
{
  size_t i;
  ps_reactant_container_t *current_rc;

  // Loop through alloc_count because containers may have been
  // initialized previously, but never reused.
  for (i = 0, current_rc = dynps_reactant_containers->ps_reactant_containers.a;
       i < dynps_reactant_containers->alloc_count;
       i++, current_rc++)
  {
    ps_reactant_container_destroy(current_rc);
  }
  free(dynps_reactant_containers->ps_reactant_containers.a);
  free(dynps_reactant_containers);
}

void
dynps_reactant_containers_output
  (const output_type ptype, const int level,
   dynps_reactant_containers_t *dynps_reactant_containers)
{
  int i;
  ps_reactant_container_t *current_ps_reactant_container;

  for (i = 0, current_ps_reactant_container =
              dynps_reactant_containers->ps_reactant_containers.a;
       i < dynps_reactant_containers->count;
       i++, current_ps_reactant_container++)
  {
    CUSTOM_PRINTF(ptype, level, " %d(%s){",
                 current_ps_reactant_container->reaction_definition->index,
                 current_ps_reactant_container->reaction_definition->name);
    RULE_AGGREGATE_OUTPUT(ptype, level,
                          current_ps_reactant_container->aggregate);
    CUSTOM_PRINTF(ptype, level, "}[%d]",
                 current_ps_reactant_container->reactant_index);
  }

  CUSTOM_PRINTF(ptype, level, "\n");
}

