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

#ifndef DYNPS_REACTANT_CONTAINERS_H
#define DYNPS_REACTANT_CONTAINERS_H

typedef struct dynps_reactant_containers_struct dynps_reactant_containers_t;
typedef union ps_reactant_containerarray_union ps_reactant_containerarray_t;

#include <stdlib.h>

#include "nauty.h"
#include "ps_reactant_container.h"
#include "reaction_definition.h"
#include "output.h"

#define DYNPS_REACTANT_CONTAINERS_OUTPUT(ptype, level, containers) \
   CHECKDEBUG(dynps_reactant_containers_output, ptype, level, containers)

union ps_reactant_containerarray_union
{
  ps_reactant_container_t *a;
  void *v;
};

struct dynps_reactant_containers_struct
{
  // Pointer to allocated array of ps_reactant_containers.
  ps_reactant_containerarray_t ps_reactant_containers;

  // Number of elements currently in use.
  size_t count;

  // Number of elements currently allocated.
  size_t alloc_count;
};

boolean
dynps_reactant_containers_ps_reactant_container_add
  (reaction_definition_t *reaction_definition, const int reactant_index,
   const int primary_particle, const size_t aggregate_size,
   const size_t full_aggregate_size, const int *aggregate,
   dynps_reactant_containers_t *dynps_reactant_containers);
ps_reactant_container_t *dynps_reactant_containers_ps_reactant_container_alloc
                        (dynps_reactant_containers_t *dynps_reactant_containers);
dynps_reactant_containers_t *dynps_reactant_containers_create
                            (const size_t init_count);
void dynps_reactant_containers_destroy(dynps_reactant_containers_t *da);
void dynps_reactant_containers_output
       (const output_type ptype, const int level,
        dynps_reactant_containers_t *dynps_reactant_containers);

#endif
