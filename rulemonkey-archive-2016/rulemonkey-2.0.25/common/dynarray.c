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


#include <limits.h>
#include <stdio.h>
#include "dynarray.h"
#include <assert.h>
#include "output.h"

void suppress_warnings25() {SUPPRESS_WARNINGS;}

boolean
dynarray_create(size_t count, size_t *palloc_count,
                size_t size, void **pdata)
{
  boolean status;

  status = dynarray_realloc(0, count, size, pdata);
  if (FALSE == status)
  {
    error_printf("Unable to create array, count = %ld, size = %ld.\n",
                 count, size);
    return FALSE;
  }

  *palloc_count = count;

  return TRUE;
}

boolean
dynarray_create2(size_t count, size_t *palloc_count,
                 size_t size1, void **pdata1,
                 size_t size2, void **pdata2)
{
  int status;

  status = dynarray_realloc(0, count, size1, pdata1);
  if (FALSE == status)
  {
    error_printf("Unable to create array, count = %ld, size = %ld.\n",
                 count, size1);
    return FALSE;
  }

  status = dynarray_realloc(0, count, size2, pdata2);
  if (FALSE == status)
  {
    free(*pdata2);
    error_printf("Unable to create array, count = %ld, size = %ld.\n",
                 count, size2);
    return FALSE;
  }

  *palloc_count = count;

  return TRUE;
}

boolean
dynarray_create3(size_t count, size_t *palloc_count,
                 size_t size1, void **pdata1,
                 size_t size2, void **pdata2,
                 size_t size3, void **pdata3)
{
  int status;

  status = dynarray_realloc(0, count, size1, pdata1);
  if (FALSE == status)
  {
    error_printf("Unable to create array, count = %ld, size = %ld.\n",
                 count, size1);
    return FALSE;
  }

  status = dynarray_realloc(0, count, size2, pdata2);
  if (FALSE == status)
  {
    free(*pdata2);
    error_printf("Unable to create array, count = %ld, size = %ld.\n",
                 count, size2);
    return FALSE;
  }

  status = dynarray_realloc(0, count, size3, pdata3);
  if (FALSE == status)
  {
    free(*pdata3);
    error_printf("Unable to create array, count = %ld, size = %ld.\n",
                 count, size3);
    return FALSE;
  }

  *palloc_count = count;

  return TRUE;
}

size_t
dynarray_check_size(size_t count, size_t alloc_count)
{
  size_t new_count;
  size_t new_alloc_count;

  assert(0 != count);

  if (count <= 0)
  {
    // Need to allocate something.
    new_count = 1;
  }
  else
  {
    new_count = count;
  }

  if (alloc_count <= 0)
  {
    new_alloc_count = 1;
  }
  else
  {
    new_alloc_count = alloc_count;
  }

  if (alloc_count >= new_count)
  {
    // Enough memory already available.
    return 0;
  }

#if AGGRESSIVE_CHECK
  if (new_alloc_count < new_count)
  {
    new_alloc_count = new_count;
  }
#else
  do
  {
#ifndef DS_NO_DOUBLECHECK
    size_t new_new_alloc_count;

    // Allocate more memory.
    new_new_alloc_count = new_alloc_count * 2;
    if ((new_new_alloc_count / 2) != new_alloc_count)
    {
      // Overflow.
     error_printf("Overflow in dynarray_realloc: count = %ld\n", count);
     return 0;
    }

    new_alloc_count = new_new_alloc_count;
#else
    new_alloc_count *= 2;
#endif
  } while (new_alloc_count < new_count);
#endif

  return new_alloc_count;
}

boolean
dynarray_resize(size_t count, size_t *palloc_count, size_t size, void **pdata)
{
  size_t new_alloc_count;
  boolean status;

  assert(0 != count);

  new_alloc_count = dynarray_check_size(count, *palloc_count);
  if (0 == new_alloc_count)
  {
    // Enought memory already available
    return TRUE;
  }

  status = dynarray_realloc(*palloc_count, new_alloc_count, size, pdata);
  if (FALSE == status)
  {
    error_printf("dynarray_resize: Unable to realloc; count = %ld, "
                 "new_alloc_count = %ld, size = %ld\n", count, new_alloc_count,
                 size);
    return FALSE;
  }

  *palloc_count = new_alloc_count;

  return TRUE;
}

boolean
dynarray_resize2(size_t count, size_t *palloc_count,
                 size_t size1, void **pdata1,
                 size_t size2, void **pdata2)
{
  size_t new_alloc_count;
  boolean status;

  new_alloc_count = dynarray_check_size(count, *palloc_count);
  if (0 == new_alloc_count)
  {
    // Enought memory already available
    return TRUE;
  }

  status = dynarray_realloc(*palloc_count, new_alloc_count, size1, pdata1);
  if (FALSE == status)
  {
    error_printf("dynarray_resize2: Unable to realloc(1); count = %ld, "
                 "new_alloc_count = %ld, size = %ld\n", count, new_alloc_count,
                 size1);
    return FALSE;
  }

  status = dynarray_realloc(*palloc_count, new_alloc_count, size2, pdata2);
  if (FALSE == status)
  {
    error_printf("dynarray_resize2: Unable to realloc(2); count = %ld, "
                 "new_alloc_count = %ld, size = %ld\n", count, new_alloc_count,
                 size2);
    dynarray_realloc(new_alloc_count, *palloc_count, size1, pdata1);
    return FALSE;
  }

  *palloc_count = new_alloc_count;

  return status;
}

boolean
dynarray_resize3(size_t count, size_t *palloc_count,
                 size_t size1, void **pdata1,
                 size_t size2, void **pdata2,
                 size_t size3, void **pdata3)
{
  size_t new_alloc_count;
  boolean status;

  new_alloc_count = dynarray_check_size(count, *palloc_count);
  if (0 == new_alloc_count)
  {
    // Enought memory already available
    return TRUE;
  }

  status = dynarray_realloc(*palloc_count, new_alloc_count, size1, pdata1);
  if (FALSE == status)
  {
    error_printf("dynarray_resize2: Unable to realloc(1); count = %ld, "
                 "new_alloc_count = %ld, size = %ld\n", count, new_alloc_count,
                 size1);
    return FALSE;
  }

  status = dynarray_realloc(*palloc_count, new_alloc_count, size2, pdata2);
  if (FALSE == status)
  {
    error_printf("dynarray_resize2: Unable to realloc(2); count = %ld, "
                 "new_alloc_count = %ld, size = %ld\n", count, new_alloc_count,
                 size2);
    dynarray_realloc(new_alloc_count, *palloc_count, size1, pdata1);
    return FALSE;
  }

  status = dynarray_realloc(*palloc_count, new_alloc_count, size3, pdata3);
  if (FALSE == status)
  {
    error_printf("dynarray_resize3: Unable to realloc(3); count = %ld, "
                 "new_alloc_count = %ld, size = %ld\n", count, new_alloc_count,
                 size3);
    dynarray_realloc(new_alloc_count, *palloc_count, size1, pdata1);
    return FALSE;
  }

  *palloc_count = new_alloc_count;

  return status;
}

boolean
dynarray_realloc(size_t old_alloc_count, size_t new_alloc_count,
                 size_t size, void **pdata)
{
  void *new_data;
  size_t old_size;
  size_t alloc_size;

  if (0 == new_alloc_count)
  {
    // Make sure we allocate something.
    new_alloc_count = 1;
  }

  old_size = old_alloc_count * size;

  alloc_size = new_alloc_count * size;
#ifndef DS_NO_DOUBLECHECK
  if ((alloc_size / new_alloc_count) != size)
  {
    // Overflow.
    error_printf("Overflow in dynarray_realloc: new_alloc_count = %ld, "
                 "size = %ld.\n", new_alloc_count, size);
    return 0;
  }
#endif

  new_data = realloc(*pdata, alloc_size);
  if (NULL == new_data)
  {
    error_printf("Unable to reallocate array, count = %ld, size = %ld.\n",
                 new_alloc_count, size);
    return 0;
  }

  if (new_alloc_count > old_alloc_count)
    {
      // Clear new memory that was allocated. */
      memset(new_data + old_size, 0, alloc_size - old_size);
    }

  *pdata = new_data;

  return 1;
}

