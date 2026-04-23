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

/*
** ps_specie.h
**
** This library provides everything needed to work with ps_species.
*/
#ifndef PS_SPECIE_H
#define PS_SPECIE_H

typedef struct ps_specie_struct ps_specie_t;

#include "dynps_reactant_containers.h"
#include "dynps_reaction_containers.h"
#include "reaction_definition.h"
#include "world.h"
#include "ps_particle.h"
#include "dynps_species.h"
#include "dynspecies.h"
#include "dynrule_aggregates.h"
#include "ps_particle.h"
#include "output.h"

#define PS_SPECIE_OUTPUT(ptype, level, ps_specie, specie) \
   CHECKDEBUG(ps_specie_output, ptype, level, ps_specie, specie)

struct ps_specie_struct
{
  // Array of reactant container lists (one list per particle) that
  // represent specie in a given reaction so that it is easy to tell
  // if two particles can react.
  dynps_reactant_containers_t **dynps_reactant_containers_array;

  // When we reuse this structure, we need to know the allocated array size.
  size_t n;
};
boolean ps_specie_setup(cgraph_t *cg, ps_specie_t *new_ps_specie);
void ps_specie_reset(ps_specie_t *ps_specie);
/* ps_specie_update() not needed. */
boolean ps_specie_check_all_reactants(ps_specie_t *ps_specie, specie_t *specie,
                                      world_t *world);
void ps_specie_output(const output_type ptype, const int level,
                      ps_specie_t *ps_specie, specie_t *specie);
ps_particle_t *ps_dynspecies_setup_particle(dynps_species_t *dynps_species,
                                     dynspecies_t *dynspecies);
boolean ps_specie_attempt_reaction
  (ps_particle_t *ps_particle1, ps_particle_t *ps_particle2,
   ps_particle_t *ps_particles, world_t *world, dynps_species_t *dynps_species,
   dynps_reaction_containers_t *dynps_reaction_containers);

#endif
