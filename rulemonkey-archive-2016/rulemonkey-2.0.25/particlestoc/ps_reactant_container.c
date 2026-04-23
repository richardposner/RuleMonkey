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
#include "ps_reactant_container.h"
#include "constants.h"

void suppress_warnings5() {SUPPRESS_WARNINGS;}

boolean
ps_reactant_container_setup
  (ps_reactant_container_t *ps_reactant_container,
   reaction_definition_t *reaction_definition, const int reactant_index,
   const int primary_particle, const size_t aggregate_size,
   const size_t full_aggregate_size, const int *aggregate)
{
  // Set structure members.
  assert(NULL != reaction_definition);
  ps_reactant_container->reaction_definition = reaction_definition;
  ps_reactant_container->reactant_index = reactant_index;
  ps_reactant_container->primary_particle = primary_particle;
  ps_reactant_container->full_aggregate_size = full_aggregate_size;

  if (ps_reactant_container->aggregate != NULL)
  {
    // We are reusing this object.
    ps_reactant_container->aggregate =
      realloc(ps_reactant_container->aggregate,aggregate_size);
  }
  else
  {
    // We are creating this object for the first time.
    ps_reactant_container->aggregate = malloc(aggregate_size);
  }
  if (NULL == ps_reactant_container->aggregate)
  {
    return FALSE;
  }
  memcpy(ps_reactant_container->aggregate, aggregate, aggregate_size);

  return TRUE;
}

void
ps_reactant_container_destroy(ps_reactant_container_t *reactant_container)
{
  if (reactant_container->aggregate != NULL)
  {
    free(reactant_container->aggregate);
    reactant_container->aggregate = NULL;
  }
}

