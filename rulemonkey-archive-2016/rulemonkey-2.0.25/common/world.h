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

#ifndef WORLD_H
#define WORLD_H

typedef struct world_struct world_t;

#define DS_FLAG_BASE 1
#define DS_FLAG_ENABLE_REJECTION (DS_FLAG_BASE << 0)
#define DS_FLAG_IGNORE_RINGS (DS_FLAG_BASE << 1)

#include "llist.h"
#include "specie.h"
#include "reaction_definition.h"
#include "dynrules.h"
#include "dynspecies.h"
#include "dyncolors.h"
#include "report_definition.h"
#include "output.h"

struct world_struct
{
  // Controls different features.
  int flags;

  // Head of linked list of reaction definitions.
  reaction_definition_t *reaction_definition_list;

  // Dynamic array or rules.
  dynrules_t *dynrules;

  // Dynamic array of string to integer mapping for colors.
  dyncolors_t *dyncolors;

  // Dynamic array of species.
  dynspecies_t *dynspecies;

  // Rules used to compose reports.
  report_definition_t *report_definition_list;

  // New report definitions always added to end of list.
  report_definition_t *report_definition_last;

  // Used when REACTION_FLAG_DIFFUSION is set in any reactions.
  double diffusion_prob;
  double diffusion_denominator;

  double current_gen;
  double current_time;
  double reaction_count;

  // Type double so that it can hold a larger number of combinations.
  double total_combinations;
  double total_weight;
};

world_t *world_create();
boolean world_complete(world_t *world);
void world_output_pop(const output_type ptype, const int level,
                      world_t *world, double log_time);

#endif
