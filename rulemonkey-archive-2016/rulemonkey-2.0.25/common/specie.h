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
** specie.h
**
** This library provides everything needed to work with species.
*/
#ifndef SPECIE_H
#define SPECIE_H

typedef struct specie_struct specie_t;

#include "cgraph.h"
#include "llist.h"
#include "world.h"
#include "dyncolors.h"
#include "dynrule_aggregates.h"
#include "output.h"
#include "dynreport_definition_containers.h"

#define SPECIE_OUTPUT(ptype, level, specie, dyncolors) \
   CHECKDEBUG(specie_output, ptype, level, specie, dyncolors)

#define SPECIE_OUTPUT_NET(ptype, level, specie, dyncolors) \
   CHECKDEBUG(specie_output_net, ptype, level, specie, dyncolors)

struct specie_struct
{
  // For reporting purposes.
  int id;

  // Optional name to identify specie.
  char *name;

  // Graph that represents specie.
  cgraph_t *cg;

  // Current population of specie.
  double pop;

  // List of report definitions that have matched this specie.
  dynreport_definition_containers_t *dynreport_definition_containers;

  // If true, don't mark for deletion so new members can be made at runtime.
  boolean preserve;

  // If false, a reaction doesn't change the population.
  boolean update;
};

boolean specie_setup(int id, cgraph_t *cg, const char *name,
                     const double pop, boolean update, boolean preserve,
                     specie_t *new_specie);
void specie_update(specie_t *specie, const double pop_delta, world_t *world);
int specie_check(specie_t *specie, dynrule_aggregates_t *dynrule_aggregates,
                 reaction_definition_t *reaction_definition,
                 int reactant_index, dynrules_t *dynrules,
                 dyncolors_t *dyncolors, int *pcount);
void specie_output(const output_type ptype, const int level,
                   specie_t *specie, dyncolors_t *dyncolors);
void specie_output_net(const output_type ptype, const int level,
                       specie_t *specie, dyncolors_t *dyncolors);

#endif
