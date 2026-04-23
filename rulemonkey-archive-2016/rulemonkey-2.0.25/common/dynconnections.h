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

#ifndef DYNCONNECTIONS_H
#define DYNCONNECTIONS_H

typedef struct dynconnections_struct dynconnections_t;

#define DYNCONNECTIONS_OUTPUT(ptype, level, ...) \
   CHECKDEBUG(dynconnections_output, ptype, level, __VA_ARGS__)

#include <stdlib.h>
#include "constants.h"
#include "cgraph.h"
#include "dynarray.h"
#include "dynints.h"
#include "dynrules.h"
#include "dynconnections.h"
#include "dyncolors.h"
#include "output.h"

struct dynconnections_struct
{
  // Pointer to allocated array of strings.
  stringarray_t names;

  // Seperate species may share connections with same name.
  // Connections with same name but different specie are still seperate.
  intarray_t specie_rule_indexes;

  // Pointer to allocated array of pointers to dynints.
  // Each dynints object contains a list of connection endpoints.
  dynintsarray_t dynints;

  // Number of elements currently in use.
  size_t count;

  // Number of elements currently allocated.
  size_t alloc_count;
};

size_t dynconnections_name_update(const char *name, const double value,
                              dynconnections_t *dynconnections);
boolean
dynconnections_connection_add(const int specie_rule_index, const char *name,
                              const int node,
                              dynconnections_t *dynconnections);
void dynconnections_output(const output_type ptype, const int level,
                           dynconnections_t *dynconnections);
boolean dynconnections_compare_rules(size_t lhs_rule_index,
                                     size_t rhs_rule_index, size_t rhs_index,
                                     dynconnections_t *dynconnections);
void dynconnections_apply_to_graph(cgraph_t *cg,
                                   dynconnections_t *dynconnections);
boolean dynconnections_apply_to_rules(const int lhs_count,
                                      const int rhs_count,
                                      int strict,
                                      const int *specie_rule_indexes,
                                      dyncolors_t *dyncolors,
                                      dynrules_t *dynrules,
                                      dynconnections_t *dynconnections);
boolean dynconnections_apply_map(const int lhs_count, const int rhs_count,
                                 const int *specie_rule_indexes,
                                 dynrules_t *dynrules, dynmaps_t *dynmaps,
                                 dynconnections_t *dynconnections);
size_t dynconnections_connection_alloc(dynconnections_t *dynconnections);
size_t dynconnections_get_index(const int specie_rule_index, const char *name,
                                dynconnections_t *dynconnections);
dynconnections_t * dynconnections_create(const size_t input_init_count);
void dynconnections_reset(dynconnections_t *dynconnections);
void dynconnections_destroy(dynconnections_t *dynconnections);

#endif
