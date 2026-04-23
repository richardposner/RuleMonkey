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


#ifndef PS_PARTICLE_H
#define PS_PARTICLE_H

typedef struct ps_particle_struct ps_particle_t;

#include "nauty.h"
#include "cgraph.h"
#include "dynps_species.h"
#include "output.h"

#define PS_PARTICLE_OUTPUT_ALL(ptype, level, ps_particles, \
                               ps_particle_count, dynps_species) \
   CHECKDEBUG(ps_particle_output_all, ptype, level, ps_particles, \
              ps_particle_count, dynps_species)

#define PS_PARTICLE_OUTPUT_SPECIE(ptype, level, ps_particle) \
   CHECKDEBUG(ps_particle_output_specie, ptype, level, ps_particle)

#define PS_PARTICLE_OUTPUT(ptype, level, ps_particle,\
                           dynps_species, ps_particles) \
   CHECKDEBUG(ps_particle_output, ptype, level, ps_particle, \
              dynps_species, ps_particles)

struct ps_particle_struct {
  // Specie that particle is currently part of.
  size_t specie_index;

  // Node index in specie.
  int specie_particle_index;

  // Keep direct link to save time later.
  dynps_reactant_containers_t *dynps_reactant_containers;

  // Circular linked list to associate all particles in current specie.
  ps_particle_t *next;
};


void
ps_particle_setup(const size_t specie_index, const int specie_particle_index,
                  dynps_reactant_containers_t *dynps_reactant_containers,
                  ps_particle_t *next, ps_particle_t *new_ps_particle);

void
ps_particle_output_all(const output_type ptype, const int level,
                       ps_particle_t *ps_particles, int ps_particle_count,
                       dynps_species_t *dynps_species);

void
ps_particle_output(const output_type ptype, const int level,
                   ps_particle_t *ps_particle, dynps_species_t *dynps_species,
                   ps_particle_t *ps_particles);

void
ps_particle_output_specie(const output_type ptype, const int level,
                          ps_particle_t *ps_particle);

boolean
ps_particle_check_same_entity(ps_particle_t *particle1,
                               ps_particle_t *particle2);

boolean
ps_particle_cgraph_append(ps_particle_t **pfull_particle,
                          cgraph_t *full_cgraph,
                          int *full_aggregate,
                          ps_particle_t *particle1,
                          cgraph_t *cgraph1,
                          int *aggregate1,
                          cgraph_t *final_cgraph);

boolean
ps_particle_add_new_species(ps_particle_t *beginning_particle,
                            cgraph_t *full_cgraph, world_t *world,
                            dynps_species_t *dynps_species);



#endif
