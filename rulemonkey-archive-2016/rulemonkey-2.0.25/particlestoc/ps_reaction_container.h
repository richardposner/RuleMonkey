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


#ifndef PS_REACTION_CONTAINER_H
#define PS_REACTION_CONTAINER_H

typedef struct ps_reaction_container_struct ps_reaction_container_t;

#include "reaction_definition.h"
#include "world.h"
#include "dynps_species.h"
#include "ps_particle.h"
#include "output.h"

#define PS_REACTION_CONTAINER_OUTPUT(ptype, level, ...) \
   CHECKDEBUG(ps_reaction_container_output, ptype, level, __VA_ARGS__)

struct ps_reaction_container_struct
{
  // Reaction to perform.
  reaction_definition_t *reaction_definition;

  // Which reactant each is for.
  int reactant_index1;
  int reactant_index2;

  // Aggregates for each reactant.
  int *aggregate1;
  int *aggregate2;

  // Size of full_aggregate when two partial aggregates are merged.
  size_t full_aggregate_size;
};

void
ps_reaction_container_setup(const size_t full_aggregate_size,
                            reaction_definition_t *reaction_definition,
                            const int reactant_index1, int *aggregate1,
                            const int reactant_index2, int *aggregate2,
                            ps_reaction_container_t *ps_reaction_container);

boolean
ps_reaction_container_perform_reaction
  (ps_reaction_container_t *ps_reaction_container,
   ps_particle_t *particle1, ps_particle_t *particle2,
   world_t *world, dynps_species_t *dynps_species);

void
ps_reaction_container_output(const output_type ptype, const int level,
                             ps_reaction_container_t *ps_reaction_container);

#endif
