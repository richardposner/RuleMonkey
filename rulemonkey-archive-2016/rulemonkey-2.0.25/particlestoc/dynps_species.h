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

#ifndef PS_DYNSPECIES_H
#define PS_DYNSPECIES_H

typedef struct dynps_species_struct dynps_species_t;
typedef union ps_speciearray_union ps_speciearray_t;

#include <stdlib.h>
#include "constants.h"
#include "specie.h"
#include "dyncolors.h"
#include "world.h"
#include "ps_specie.h"
#include "ps_particle.h"
#include "dynrule_aggregates.h"
#include "output.h"

#define DYNPS_SPECIES_OUTPUT(ptype, level, dynps_species, \
                             dynspecies, dyncolors) \
   CHECKDEBUG(dynps_species_output, ptype, level, dynps_species, \
              dynspecies, dyncolors)

union ps_speciearray_union
{
  ps_specie_t *a;
  void *v;
};

struct dynps_species_struct
{
  // List of allocated but currently unused entries.
  ilistarray_t available;
  
  // Pointer to allocated array of ps_species.
  ps_speciearray_t ps_species;

  // Number of elements currently in use.
  size_t count;

  // Number of elements currently allocated.
  size_t alloc_count;
};

ps_particle_t *dynps_species_setup_particles
                 (dynspecies_t *dynspecies, size_t *pps_particle_count,
                  dynps_species_t *dynps_species, dyncolors_t *dyncolors);
boolean dynps_species_specie_add(cgraph_t *cg, const char *name,
                                 const double pop_delta, boolean update,
                                 boolean last, boolean match, int *pspecie_id,
                                 boolean *pnew,
                                 dynspecies_t *dynspecies, world_t *world,
                                 dynps_species_t *dynps_species);
boolean dynps_species_specie_remove
          (specie_t *specie, dynps_species_t *dynps_species);
boolean dynps_species_specie_create_all(dynspecies_t *dynspecies,
                                        dynps_species_t *dynps_species);
boolean dynps_species_check_all_reactants
          (dynps_species_t *dynps_species,
           dynspecies_t *dynspecies,
           world_t *world);
void dynps_species_output(const output_type ptype, const int level,
                          dynps_species_t *dynps_species,
                          dynspecies_t *dynspecies, dyncolors_t *dyncolors);
boolean dynps_species_ps_specie_alloc(dynps_species_t *dynspecies);
dynps_species_t *dynps_species_create(const size_t init_count);
void dynps_species_destroy(dynps_species_t *dynps_species);

#endif
