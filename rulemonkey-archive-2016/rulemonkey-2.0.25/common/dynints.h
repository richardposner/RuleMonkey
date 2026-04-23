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

#ifndef DYNINTS_H
#define DYNINTS_H

typedef struct dynints_struct dynints_t;
typedef union dynintsarray_union dynintsarray_t;

#include <stdlib.h>

#include "nauty.h"
#include "dynarray.h"

union dynintsarray_union
{
  dynints_t **a;
  void *v;
};

struct dynints_struct
{
  // Pointer to allocated array of ints.
  intarray_t ints;

  // Number of elements currently in use.
  size_t count;

  // Number of elements currently allocated.
  size_t alloc_count;
};

boolean dynints_int_add(const int value, dynints_t *dynints);
boolean dynints_copy(dynints_t *src, dynints_t *dest);
boolean dynints_int_exists(const int value, dynints_t *dynints);
int dynints_int_find(const int value, dynints_t *dynints);
boolean dynints_compare(dynints_t *a_states, dynints_t *b_states);
boolean dynints_compare_cond(dynints_t *a_states, dynints_t *b_states);
boolean dynints_int_remove(const int value, dynints_t *dynints);
int *dynints_int_alloc(dynints_t *dynints);
dynints_t *dynints_create(const size_t init_count);
boolean dynints_resize(dynints_t *dynints, const size_t count);
void dynints_destroy(dynints_t *da);

#endif
