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

#ifndef DYNSYMBOLS_H
#define DYNSYMBOLS_H

typedef struct dynsymbols_struct dynsymbols_t;

#include <stdlib.h>
#include "constants.h"
#include "dynarray.h"
#include "output.h"

struct dynsymbols_struct
{
  // Pointer to allocated array of strings.
  stringarray_t names;

  // Pointer to allocated array of values.
  doublearray_t values;

  // Number of elements currently in use.
  size_t count;

  // Number of elements currently allocated.
  size_t alloc_count;
};

size_t dynsymbols_name_update(const char *name, const double value,
                              dynsymbols_t *dynsymbols);
void dynsymbols_output(const output_type ptype, const int level,
                       dynsymbols_t *dynsymbols, const char *prefix,
                       const char *postfix, boolean value,
                       boolean prompt);
size_t dynsymbols_symbol_alloc(dynsymbols_t *dynsymbols);
size_t dynsymbols_get_index(const char *name, dynsymbols_t *dynsymbols);
size_t dynsymbols_get_value(const char *name, double *value,
                            dynsymbols_t *dynsymbols);
dynsymbols_t * dynsymbols_create(const size_t input_init_count);
void dynsymbols_destroy(dynsymbols_t *dynsymbols);

#endif
