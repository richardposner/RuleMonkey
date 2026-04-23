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

#ifndef DYNPS_REACTION_CONTAINERS_H
#define DYNPS_REACTION_CONTAINERS_H

typedef struct dynps_reaction_containers_struct dynps_reaction_containers_t;

#include <stdlib.h>
#include "ps_reaction_container.h"
#include "dynps_reaction_containers.h"
#include "ps_particle.h"
#include "output.h"

typedef union ps_reaction_containerarray_union ps_reaction_containerarray_t;

#define DYNPS_REACTION_CONTAINERS_OUTPUT(ptype, level, containers) \
   CHECKDEBUG(dynps_reaction_containers_output, ptype, level, containers)

union ps_reaction_containerarray_union
{
  ps_reaction_container_t *a;
  void *v;
};

struct dynps_reaction_containers_struct
{
  // Pointer to allocated array of reaction_containers.
  ps_reaction_containerarray_t ps_reaction_containers;

  // Pointer to allocated array of cumulative probabilities.
  doublearray_t cumulative_probabilities;

  // Sum of all probabilities.
  double prob_sum;

  // Number of elements currently in use.
  size_t count;

  // Number of elements currently allocated.
  size_t alloc_count;
};

void
dynps_reaction_containers_reset
  (dynps_reaction_containers_t *dynps_reaction_containers);

void
dynps_reaction_containers_output
  (const output_type ptype, const int level,
   dynps_reaction_containers_t *dynps_reaction_containers);

boolean dynps_reaction_containers_reaction_container_add
  (dynps_reaction_containers_t *dynps_reaction_containers,
   reaction_definition_t *reaction_definition,
   const int reactant_index1, int *aggregate1,
   const int reactant_index2, int *aggregate2,
   const int full_aggregate_size);

boolean dynps_reaction_containers_attempt_reaction
  (dynps_reaction_containers_t *dynps_reaction_containers,
   ps_particle_t *particle1, ps_particle_t *particle2,
   world_t *world, dynps_species_t *dynps_species);

dynps_reaction_containers_t *dynps_reaction_containers_create
  (const size_t init_count);

void dynps_reaction_containers_destroy
  (dynps_reaction_containers_t *dynps_reaction_containers);
#endif

