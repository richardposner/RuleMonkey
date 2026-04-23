/*
** Copyright 2008 Northern Arizona Univeristy and
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
** For each specie, we need to keep track of each
** reaction definition that specie participates in so that
** the cached populations can be updated or removed.
*/
#ifndef RS_SPECIE_H
#define RS_SPECIE_H

typedef struct rs_specie_struct rs_specie_t;
typedef struct rs_specie_output_d_struct rs_specie_output_d_t;

#include "specie.h"
#include "world.h"
#include "rs_reaction_definition.h"

#define RS_SPECIE_FSA_OUTPUT(ptype, level, rs_specie_fsa, world) \
   CHECKDEBUG(rs_specie_fsa_output, ptype, level, rs_specie_fsa, world)
#define RS_SPECIE_OUTPUT(ptype, level, rs_specie, world) \
   CHECKDEBUG(rs_specie_fsa_output, ptype, level, rs_specie, world)

#define RS_SPECIES_CREATE() \
   FSA_CREATE(sizeof(rs_specie_t), 8, rs_specie_init, rs_specie_clean)

struct rs_specie_struct
{
  int specie_id;

  // List of pointers to rs specie containers.
  // These are used to update specie populations in each reaction definition.
  // The fsa needs to be deleted when cleaned up, but do not cleanup contents
  // because the contents are allocated in rs_reaction_definition.
  fsa_t *rs_specie_container_pointers_fsa;
};

struct rs_specie_output_d_struct
{
  output_type ptype;
  int level;
  world_t *world;
};

int rs_specie_init(void *d, void *elem);

int rs_specie_clean(void *d, void *elem);

double rs_specie_update
         (rs_specie_t *rs_specie, specie_t *specie, double pop_delta,
          world_t *world, double *ptotal_prob);

int rs_specie_add_new_species
  (fsa_t *rs_specie_fsa, cgraph_t *full_cgraph,
   world_t *world, double *pdelta_prob,
   rs_reaction_definition_t *rs_reaction_defintions);

void rs_specie_remove(rs_specie_t *rs_specie, fsa_t *rs_specie_fsa);

int rs_specie_check_all_reactants
  (rs_specie_t *rs_specie,
   rs_reaction_definition_t *rs_reaction_definitions, specie_t *specie,
   world_t *world, double *ptotal_prob);

int rs_specie_output_cb(void *d, void *elem);
void rs_specie_output(const output_type ptype, const int level,
                      rs_specie_t *rs_specie, world_t *world);

void rs_specie_fsa_output(const output_type ptype, const int level,
                          fsa_t *rs_specie_fsa, world_t *world);

#endif
