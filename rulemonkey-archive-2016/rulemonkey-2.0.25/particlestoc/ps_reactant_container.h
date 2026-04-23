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


#ifndef PS_REACTANT_CONTAINER_H
#define PS_REACTANT_CONTAINER_H

typedef struct ps_reactant_container_struct ps_reactant_container_t;

#include "dynps_reactant_containers.h"
#include "reaction_definition.h"

#include "reaction_definition.h"

struct ps_reactant_container_struct
{
  /* Used to determine if two given particles can interact. */

  // Pointer to reaction definition that reactant participates in.
  reaction_definition_t *reaction_definition;

  // Number of reactants current reaction requires.

  // Index of particular reactant that represents specie.
  int reactant_index;

  // Primary particle used to corrolate split or self bind pairs.
  int primary_particle;

  /* Used to actually perform reaction. */

  // Aggregate containing all the actions to perform.
  int *aggregate;

  // Required size of full aggregate.
  int full_aggregate_size;
};

boolean
ps_reactant_container_setup
  (ps_reactant_container_t *ps_reactant_container,
   reaction_definition_t *reaction_definition, const int reactant_index,
   const int primary_particle, const size_t aggregate_size,
   const size_t full_aggregate_size, const int *aggregate);
void ps_reactant_container_destroy(ps_reactant_container_t *reactant_container);

#endif
