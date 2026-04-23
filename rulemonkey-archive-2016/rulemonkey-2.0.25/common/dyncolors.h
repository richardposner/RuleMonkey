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

#ifndef DYNCOLORS_H
#define DYNCOLORS_H

typedef struct dyncolors_struct dyncolors_t;

#include <stdlib.h>
#include "constants.h"
#include "color.h"
#include "dynspecies.h"
#include "world.h"
#include "output.h"

typedef union colorarray_union colorarray_t;

#define DYNCOLORS_OUTPUT(ptype, level, ...) \
   CHECKDEBUG(dyncolors_output, ptype, level, __VA_ARGS__)

union colorarray_union
{
  color_t *a;
  void *v;
};

struct dyncolors_struct
{
  // Pointer to allocated array of colors.
  colorarray_t colors;

  // Number of elements currently in use.
  size_t count;

  // Number of elements currently allocated.
  size_t alloc_count;
};

color_t *dyncolors_color_create(const char *name, color_type_t orig_type,
                                dyncolors_t *dyncolors);
color_t * real_dyncolors_color_create(const char *name, color_type_t type,
                                      dynspecies_t *dynspecies,
                                      dyncolors_t *dyncolors, world_t *world);
color_t *dyncolors_color_create_any(const char *name, color_type_t type,
                                    dynspecies_t *dynspecies,
                                    dyncolors_t *dyncolors, world_t *world);
color_t *dyncolors_color_alloc(dyncolors_t *dyncolors);
void dyncolors_output(const output_type ptype, const int level,
                      dyncolors_t *dyncolors, const char *prefix,
                      const char *postfix, boolean value, boolean prompt);
color_t *dyncolors_color_normalized_name_lookup(const char *normalized_name,
                                                dyncolors_t *dyncolors);
color_t *dyncolors_color_name_type_lookup(const char *name, color_type_t type,
                                          dyncolors_t *dyncolors);
color_t *dyncolors_color_value_lookup(const int value,
                                      dyncolors_t *dyncolors);
color_type_t dyncolors_get_type_from_value(const int value,
                                           dyncolors_t *dyncolors);
const char * dyncolors_get_name_from_value(const int value,
                                           dyncolors_t *dyncolors);
dyncolors_t *dyncolors_create(const size_t init_count);
void dyncolors_destroy(dyncolors_t *da);

#endif
