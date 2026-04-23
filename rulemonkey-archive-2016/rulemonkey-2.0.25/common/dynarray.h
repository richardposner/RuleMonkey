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

#ifndef DYNARRAY_H
#define DYNARRAY_H

typedef union intarray_union intarray_t;
typedef union stringarray_union stringarray_t;
typedef union sizearray_union sizearray_t;
typedef union doublearray_union doublearray_t;

#include <stdlib.h>
#include "constants.h"

// To prevent type-punning warnings.
union intarray_union
{
  int *a;
  void *v;
};

union stringarray_union
{
  char **a;
  void *v;
};

union sizearray_union
{
  size_t *a;
  void *v;
};

union doublearray_union
{
  double *a;
  void *v;
};

boolean dynarray_create(size_t new_alloc_count, size_t *palloc_count,
                        size_t size, void **pdata);
boolean dynarray_create2(size_t new_alloc_count, size_t *palloc_count,
                         size_t size1, void **pdata1,
                         size_t size2, void **pdata2);
boolean dynarray_create3(size_t count, size_t *palloc_count,
                         size_t size1, void **pdata1,
                         size_t size2, void **pdata2,
                         size_t size3, void **pdata3);
boolean dynarray_resize(size_t new_alloc_count, size_t *palloc_count,
                        size_t size, void **pdata);
boolean dynarray_resize2(size_t new_alloc_count, size_t *palloc_count,
                         size_t size1, void **pdata1,
                         size_t size2, void **pdata2);
boolean dynarray_resize3(size_t count, size_t *palloc_count,
                         size_t size1, void **pdata1,
                         size_t size2, void **pdata2,
                         size_t size3, void **pdata3);
size_t dynarray_check_size(size_t count, size_t alloc_count);
boolean dynarray_realloc(size_t old_alloc_count, size_t new_alloc_count,
                         size_t size, void **pdata);

#endif
