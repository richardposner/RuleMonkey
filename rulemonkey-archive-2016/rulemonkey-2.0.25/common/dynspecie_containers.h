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

#ifndef DYNSPECIE_CONTAINERS_H
#define DYNSPECIE_CONTAINERS_H

typedef struct dynspecie_containers_struct dynspecie_containers_t;
typedef union specie_containerarray_union specie_containerarray_t;

#include "dynrules.h"
#include "rule.h"
#include "specie_container.h"
#include "dyncolors.h"
#include "world.h"

union specie_containerarray_union
{
  specie_container_t *a;
  void *v;
};

struct dynspecie_containers_struct
{
  // Pointer to allocated array of specie_containers.
  specie_containerarray_t specie_containers;

  // TRUE if all specie containers are unique
  // (specie_containers[i].count == 1 for all i).
  boolean singular;

  // Number of elements currently in use for reaction.
  // (Number of nodes matched for report_definition.)
  size_t count;

  // Number of elements currently allocated.
  size_t alloc_count;
};

boolean dynspecie_containers_specie_container_add
          (specie_t *specie, dynspecie_containers_t *dynspecie_containers,
           world_t *world);
boolean dynspecie_containers_specie_container_add_count
          (specie_t *specie, const int new_count,
           dynspecie_containers_t *dynspecie_containers, world_t *world);
void dynspecie_containers_remove_count
       (dynspecie_containers_t *dynspecie_containers, world_t *world);
specie_container_t *dynspecie_containers_specie_container_alloc
                      (dynspecie_containers_t *dynspecie_containers);
dynspecie_containers_t *dynspecie_containers_create(const size_t init_count);
void dynspecie_containers_destroy
       (dynspecie_containers_t *dynspecie_containers);

#endif
