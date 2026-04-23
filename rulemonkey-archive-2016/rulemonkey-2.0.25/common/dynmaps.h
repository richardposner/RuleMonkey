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

#ifndef DYNMAPS_H
#define DYNMAPS_H

typedef struct dynmaps_struct dynmaps_t;

#define DYNMAPS_OUTPUT(ptype, level, ...) \
   CHECKDEBUG(dynmaps_output, ptype, level, __VA_ARGS__)


#include <stdlib.h>

#include "nauty.h"
#include "dynrules.h"
#include "output.h"
#include "dynarray.h"

struct dynmaps_struct
{
  // Pointer to allocated array of map between sources and destinations.
  // sources[i] == lhs_molecule_index.
  // destinations[i] == corresponding_rhs_molecule_index.
  sizearray_t sources;
  sizearray_t destinations;

  // Number of elements currently in use.
  size_t count;

  // Number of elements currently allocated.
  size_t alloc_count;
};

boolean dynmaps_map_add(const size_t source, const size_t dest,
                        dynmaps_t *dynmaps);
void dynmaps_output(const output_type ptype, const int level,
                    dynmaps_t *dynmaps);
boolean dynmaps_apply_to_rules(dynmaps_t *dynmaps, dynrules_t *dynrules);
boolean dynmaps_dest_from_source(size_t *pdest, const size_t source, 
                                 dynmaps_t *dynmaps);
boolean dynmaps_source_from_dest(size_t *psource, const size_t dest, 
                                 dynmaps_t *dynmaps);
boolean dynmaps_dest_exists(const size_t dest, dynmaps_t *dynmaps);
boolean dynmaps_map_alloc(size_t *pnew_index, dynmaps_t *dynmaps);
dynmaps_t *dynmaps_create(const size_t init_count);
void dynmaps_destroy(dynmaps_t *da);

#endif
