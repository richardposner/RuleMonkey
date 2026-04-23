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

#ifndef DYNSPECIES_H
#define DYNSPECIES_H

typedef struct dynspecies_struct dynspecies_t;

#include <stdlib.h>
#include "constants.h"
#include "specie.h"
#include "dyncolors.h"
#include "world.h"
#include "output.h"

typedef union speciearray_union speciearray_t;

#define DYNSPECIES_OUTPUT(ptype, level, ...) \
   CHECKDEBUG(dynspecies_output, ptype, level, __VA_ARGS__)

union speciearray_union
{
  specie_t *a;
  void *v;
};

struct dynspecies_struct
{
  // Pointer to allocated array of species.
  speciearray_t species;

  // Number of elements currently in use.
  size_t count;

  // Number of elements currently allocated.
  size_t alloc_count;
};

boolean dynspecies_specie_add(cgraph_t *cg, const char *name,
                              const double pop_delta, boolean update,
                              boolean last, boolean match, boolean preserve,
                              int *pspecie_id, boolean *pnew,
                              dynspecies_t *dynspecies, world_t *world);
boolean dynspecies_specie_create(cgraph_t *cg, const char *name,
                                 const double pop, boolean update,
                                 boolean preserve,
                                 dynspecies_t *dynspecies);
void dynspecies_output(const output_type ptype, const int level,
                       dynspecies_t *dynspecies, dyncolors_t *dyncolors);
void dynspecies_output_net(const output_type ptype, const int level,
                           dynspecies_t *dynspecies, dyncolors_t *dyncolors);
void dynspecies_output_pop(const output_type ptype, const int level,
                           dynspecies_t *dynspecies, world_t *world);
boolean dynspecies_specie_alloc(dynspecies_t *dynspecies);
dynspecies_t *dynspecies_create(const size_t init_count);
void dynspecies_destroy(dynspecies_t *dynspecies);

#endif
