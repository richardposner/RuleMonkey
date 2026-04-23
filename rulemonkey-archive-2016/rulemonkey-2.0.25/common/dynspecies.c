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
#include "dynspecies.h"
#include "dynarray.h"
#include "output.h"

void suppress_warnings22() {SUPPRESS_WARNINGS;}

// Local function declarations.
static boolean dynspecies_resize(dynspecies_t *dynspecies, const size_t count);

/*
** dynspecies_specie_add
**
** Purpose: If colored graph already exists, simply add to existing
**   population.  If colored graph does not exist, create new specie
**   in population.
**
**
** Note: cg must be in a normalized form.
**
** Variables:
**   cg: Canonical & normalized colored graph representation of specie.
**   count: Number of specie to add.
**   update: Whether or not specie should be updated when reaction occurs.
**   last: If adding multiple species, FALSE for all except TRUE for last.
**   match: If FALSE, must call specie_match() later.
**   pspecie_id: Will contain new specie_id if successfyl.
**   pnew: Pointer to flag.  TRUE if new specie, otherwise will be FALSE.
**   dynspecies: Dynamic list of species.
**   world: Current state information of the simulated world.
*/
boolean
dynspecies_specie_add(cgraph_t *cg, const char *name,
                      const double pop_delta, boolean update,
                      boolean last, boolean match, boolean preserve,
                      int *pspecie_id, boolean *pnew,
                      dynspecies_t *dynspecies, world_t *world)
{
  int status;
  specie_t *new_specie;

  // New specie, add to list.
  *pnew = TRUE;
  status = dynspecies_specie_create(cg, name, pop_delta, update, preserve,
                                    dynspecies);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

  *pspecie_id = world->dynspecies->count - 1;
  new_specie = world->dynspecies->species.a + *pspecie_id;

  SPECIE_OUTPUT(OUTPUT_DEBUG, 3, new_specie, world->dyncolors);

  return status;
}

boolean
dynspecies_specie_create(cgraph_t *cg, const char *name, const double pop,
                         boolean update, boolean preserve,
                         dynspecies_t *dynspecies)
{
  boolean status;
  int index;

  // Need to create new specie.
  status = dynspecies_specie_alloc(dynspecies);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

  index = dynspecies->count - 1;
  status = specie_setup(index, cg, name, pop, update, preserve,
                        dynspecies->species.a + index);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

  return TRUE;
}

void
dynspecies_output_pop(const output_type ptype, const int level,
                      dynspecies_t *dynspecies, world_t *world)
{
  size_t i;
  specie_t *current_specie;
  size_t count;

  count = dynspecies->count;

  for (i = 0, current_specie = dynspecies->species.a;
       i < count;
       i++, current_specie++)
  {
    CUSTOM_PRINTF(ptype, level, "\t%.0f", current_specie->pop);
  }

  // End line.
  CUSTOM_PRINTF(ptype, level, "\n");
}

void
dynspecies_output(const output_type ptype, const int level,
                  dynspecies_t *dynspecies, dyncolors_t *dyncolors)
{
  size_t i;
  specie_t *current_specie;
  size_t count;

  count = dynspecies->count;
  for (i = 0, current_specie = dynspecies->species.a;
       i < count;
       i++, current_specie++)
  {
    SPECIE_OUTPUT(ptype, level, current_specie, dyncolors);
  }
}

void
dynspecies_output_net(const output_type ptype, const int level,
                      dynspecies_t *dynspecies, dyncolors_t *dyncolors)
{
  size_t i;
  specie_t *current_specie;
  size_t count;

  count = dynspecies->count;
  for (i = 0, current_specie = dynspecies->species.a;
       i < count;
       i++, current_specie++)
  {
    SPECIE_OUTPUT_NET(ptype, level, current_specie, dyncolors);
  }
}

boolean
dynspecies_specie_alloc(dynspecies_t *dynspecies)
{
  int status;

  status =
    dynspecies_resize(dynspecies, dynspecies->count + 1);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

  dynspecies->count++;

  return TRUE;
}

dynspecies_t *
dynspecies_create(const size_t init_count)
{
  dynspecies_t *new_dynspecies;
  int status;

  // All fields need to be initialized to zero anyways, so use calloc.
  new_dynspecies = calloc(1, sizeof(dynspecies_t));
  if (NULL == new_dynspecies)
  {
    error_printf("Unable to allocate new dynspecies.\n");
    return NULL;
  }

  status = dynarray_create(init_count,
                           &new_dynspecies->alloc_count, sizeof(specie_t),
                           &new_dynspecies->species.v);
  if (0 == status)
  {
    // Error message output in dynarray_create.
    free(new_dynspecies);
    return NULL;
  }

  return new_dynspecies;
}

static boolean
dynspecies_resize(dynspecies_t *dynspecies,
                  const size_t count)
{
  int status;

  status = dynarray_resize(count, &dynspecies->alloc_count,
                           sizeof(specie_t), &dynspecies->species.v);
  if (0 == status)
  {
    // Error message output in dynarray_resize.
    return FALSE;
  }

  return TRUE;
}

void
dynspecies_destroy(dynspecies_t *dynspecies)
{
  free(dynspecies->species.a);
  free(dynspecies);
}
