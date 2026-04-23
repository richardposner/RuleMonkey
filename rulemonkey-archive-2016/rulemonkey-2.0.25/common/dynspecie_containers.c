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
#include "dynspecie_containers.h"
#include "dynarray.h"
#include "constants.h"
#include "output.h"

void suppress_warnings19() {SUPPRESS_WARNINGS;}

// Local function declarations.
static int dynspecie_containers_resize
             (dynspecie_containers_t *dynspecie_containers, const size_t count);
static boolean dynspecie_containers_specie_container_create
                 (specie_t *specie, const int count,
                  dynspecie_containers_t *dynspecie_containers);

static boolean
dynspecie_containers_specie_container_create
  (specie_t *specie, const int count,
   dynspecie_containers_t *dynspecie_containers)
{
  specie_container_t *new_specie_container;

  // Need to create new specie_container.
  new_specie_container =
    dynspecie_containers_specie_container_alloc(dynspecie_containers);
  if (NULL == new_specie_container)
  {
    assert(0);
    return FALSE;
  }

  new_specie_container->specie_id = specie->id;
  new_specie_container->count = count;

  return TRUE;
}

boolean
dynspecie_containers_specie_container_add
  (specie_t *specie, dynspecie_containers_t *dynspecie_containers,
   world_t *world)
{
  int status;

  status = dynspecie_containers_specie_container_add_count
             (specie, 1, dynspecie_containers, world);

  return status;
}

boolean
dynspecie_containers_specie_container_add_count
  (specie_t *specie, const int new_count,
   dynspecie_containers_t *dynspecie_containers, world_t *world)
{
  int i;
  int count;
  boolean status;
  specie_t *current_specie;
  specie_container_t *current_specie_container;

  count = dynspecie_containers->count;
  for (i = 0, current_specie_container =
                dynspecie_containers->specie_containers.a;
       i < count;
       i++, current_specie_container++)
  {
    if (current_specie_container->specie_id == specie->id)
    {
      // Specie already present.
      dynspecie_containers->singular = FALSE;
      current_specie_container->count += new_count;

      current_specie =
        world->dynspecies->species.a + current_specie_container->specie_id;

      return TRUE;
    }
  }

  // Specie not present.
  status = dynspecie_containers_specie_container_create
             (specie, new_count, dynspecie_containers);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

  return TRUE;
}

void
dynspecie_containers_remove_count(dynspecie_containers_t *dynspecie_containers,
                                  world_t *world)
{
  int i;
  int specie_containers_count;
  specie_container_t *current_specie_container;
  specie_t *current_specie;

  // Perform last update outside of loop so that flag can be changed.
  specie_containers_count = dynspecie_containers->count;
  if (0 == specie_containers_count)
  {
    // Nothing to do.
    return;
  }
  else if ( 1 == specie_containers_count)
  {
    // Only one specie to update.
    current_specie =
      world->dynspecies->species.a +
      dynspecie_containers->specie_containers.a->specie_id;

    specie_update(current_specie,
                  -(double)dynspecie_containers->specie_containers.a->count,
                  world);

    return;
  }

  // More than one specie to update.  Update last specie seperately for speed.
  for (i = 0, current_specie_container =
                dynspecie_containers->specie_containers.a;
       i < specie_containers_count - 1;
       i++, current_specie_container++)
  {
    current_specie =
      world->dynspecies->species.a + current_specie_container->specie_id;

    if (TRUE == current_specie->update)
    {
      specie_update(current_specie,
                    -(double)current_specie_container->count, world);
    }
  }

  // Perform final update outside of loop.
  current_specie =
    world->dynspecies->species.a + current_specie_container->specie_id;

  specie_update(current_specie,
                -(double)current_specie_container->count, world);
}

specie_container_t *
dynspecie_containers_specie_container_alloc(dynspecie_containers_t *dynspecie_containers)
{
  int status;

  status = dynspecie_containers_resize(dynspecie_containers,
                                       dynspecie_containers->count + 1);
  if (FALSE == status)
  {
    assert(0);
    return NULL;
  }

  return &dynspecie_containers->specie_containers.a[dynspecie_containers->count++];
}

dynspecie_containers_t *
dynspecie_containers_create(const size_t init_count)
{
  dynspecie_containers_t *new_dynspecie_containers;
  int status;

  // All fields need to be initialized to zero anyways, so use calloc.
  new_dynspecie_containers = calloc(1, sizeof(dynspecie_containers_t));
  if (NULL == new_dynspecie_containers)
  {
    error_printf("Unable to allocate new dynspecie_containers.\n");
    return NULL;
  }

  new_dynspecie_containers->singular = TRUE;

  status = dynarray_create(init_count,
                           &new_dynspecie_containers->alloc_count,
                           sizeof(specie_container_t),
                           &new_dynspecie_containers->specie_containers.v);
  if (0 == status)
  {
    // Error message output in dynarray_create.
    free(new_dynspecie_containers);
    return NULL;
  }

  return new_dynspecie_containers;
}

static int
dynspecie_containers_resize(dynspecie_containers_t *dynspecie_containers,
                            const size_t count)
{
  int status;

  status =
    dynarray_resize(count,
                    &dynspecie_containers->alloc_count,
                    sizeof(specie_container_t),
                    &dynspecie_containers->specie_containers.v);
  if (0 == status)
  {
    // Error message output in dynarray_resize.
    return 0;
  }

  return 1;
}

void
dynspecie_containers_destroy(dynspecie_containers_t *dynspecie_containers)
{
  int i;
  specie_container_t *current_specie_container;

  for (i = 0,
        current_specie_container = dynspecie_containers->specie_containers.a;
       i <dynspecie_containers->count;
       i++, current_specie_container++)
  {
    //specie_container_destroy(current_specie_container);
  }

  free(dynspecie_containers->specie_containers.a);
  free(dynspecie_containers);
}
