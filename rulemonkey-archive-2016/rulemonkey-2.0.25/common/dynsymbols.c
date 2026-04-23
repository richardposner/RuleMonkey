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


#include <assert.h>
#include "dynsymbols.h"
#include "output.h"

void suppress_warnings23() {SUPPRESS_WARNINGS;}

// Local function declarations.
static boolean dynsymbols_resize(dynsymbols_t *dynsymbols,
                                 const size_t count);

size_t
dynsymbols_name_update(const char *name, const double value,
                       dynsymbols_t *dynsymbols)
{
  size_t i;

  i = dynsymbols_get_index(name, dynsymbols);
  if (0 != i)
  {
    // Color already exists.
    dynsymbols->values.a[i] = value;
    return i;
  }

  // Need to create new symbol.
  i = dynsymbols_symbol_alloc(dynsymbols);
  if (0 == i)
  {
    return 0;
  }

  dynsymbols->names.a[i] = strdup(name);
  if (NULL == dynsymbols->names.a[i])
  {
    return 0;
  }

  dynsymbols->values.a[i] = value;

  return i;
}

void
dynsymbols_output(const output_type ptype, const int level,
                  dynsymbols_t *dynsymbols, const char *prefix,
                  const char *postfix, boolean value,
                  boolean prompt)
{
  size_t i;
  size_t count;
  char **names;
  double *values;

  count = dynsymbols->count;
  names = dynsymbols->names.a;
  values = dynsymbols->values.a;
  // 0 is reserved to signify an error.
  for (i = 1; i < count; i++)
  {
    CUSTOM_PRINTF(ptype, level, "%s%s = %f%s",
                  prefix, names[i], values[i], postfix);
  }

  CUSTOM_PRINTF(ptype, level, "\n");
}

size_t
dynsymbols_symbol_alloc(dynsymbols_t *dynsymbols)
{
  int status;

  status = dynsymbols_resize(dynsymbols, dynsymbols->count + 1);
  if (FALSE == status)
  {
    assert(0);
    return 0;
  }

  return dynsymbols->count++;
}

size_t
dynsymbols_get_index(const char *name, dynsymbols_t *dynsymbols)
{
  size_t i;
  size_t count;
  char **names;

  if (NULL == name)
  {
    assert(0 && "dynsymbols_symbol_lookup got an empty name pointer");
    return 0;
  }

  names = dynsymbols->names.a;
  count = dynsymbols->count;
  // Index 0 is reserved for errors.
  for (i = 1; i < count; i++)
  {
    if (0 == strcasecmp(name, names[i]))
    {
      // Found a match.
      return i;
    }
  }

  // Name was not found.
  return 0;
}

size_t
dynsymbols_get_value(const char *name, double *value,
                     dynsymbols_t *dynsymbols)
{
  int i;

  i = dynsymbols_get_index(name, dynsymbols);
  *value = dynsymbols->values.a[i];

  return i;
}

dynsymbols_t *
dynsymbols_create(const size_t input_init_count)
{
  dynsymbols_t *new_dynsymbols;
  int status;
  size_t init_count;

  if (0 == input_init_count)
  {
    // Index 0 is reserved for errors, make sure it exists.
    init_count = 1;
  }
  else
  {
    init_count = input_init_count;
  }

  // All fields need to be initialized to zero anyways, so use calloc.
  new_dynsymbols = calloc(1, sizeof(dynsymbols_t));
  if (NULL == new_dynsymbols)
  {
    error_printf("Unable to allocate new dynsymbols.\n");
    return NULL;
  }

  status = dynarray_create2(init_count,
                            &new_dynsymbols->alloc_count,
                            sizeof(char *), &new_dynsymbols->names.v,
                            sizeof(double), &new_dynsymbols->values.v);
  if (0 == status)
  {
    // Error message output in dynarray_create2.
    free(new_dynsymbols);
    return NULL;
  }

  // Reserve 0 index for errors.
  new_dynsymbols->count = 1;
  new_dynsymbols->names.a[0] = "myerror";
  new_dynsymbols->values.a[0] = 0.1;

  return new_dynsymbols;
}

static boolean
dynsymbols_resize(dynsymbols_t *dynsymbols, const size_t count)
{
  int status;

  status =
    dynarray_resize2(count,
                     &dynsymbols->alloc_count,
                     sizeof(char *), &dynsymbols->names.v,
                     sizeof(double), &dynsymbols->values.v);
  if (0 == status)
  {
    // Error message output in dynarray_resize2.
    return 0;
  }

  return 1;
}

void
dynsymbols_destroy(dynsymbols_t *dynsymbols)
{
  size_t i;
  char **names;
  size_t count;

  names = dynsymbols->names.a;
  count = dynsymbols->count;
  // Index 0 is reserved for errors.
  for (i = 1; i < count; i++)
  {
    free(names[i]);
  }

  free(dynsymbols->names.a);
  free(dynsymbols->values.a);
  free(dynsymbols);
}

