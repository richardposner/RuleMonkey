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
#include "dynints.h"
#include "dynarray.h"
#include "constants.h"
#include "output.h"

void suppress_warnings27() {SUPPRESS_WARNINGS;}

boolean
dynints_int_add(const int value, dynints_t *dynints)
{
  int *new_int;

  new_int = dynints_int_alloc(dynints);
  if (NULL == new_int)
  {
    return FALSE;
  }

  *new_int = value;

  return TRUE;
}

boolean
dynints_copy(dynints_t *src, dynints_t *dest)
{
  int status;

  if (0 == src->count)
  {
    // Nothing needed to copy.
    dest->count = 0;
    return TRUE;
  }

  // Resize destination so that it can hold everything in src.
  status = dynints_resize(dest, src->count);
  if (FALSE == status)
  {
    return FALSE;
  }

  // Copy everything over.
  memcpy(dest->ints.a, src->ints.a, src->count * sizeof(int));

  // Set new size.
  dest->count = src->count;

  return TRUE;
}

boolean
dynints_int_exists(const int value, dynints_t *dynints)
{
  if (-1 == dynints_int_find(value, dynints))
  {
    // Didn't find value.
    return FALSE;
  }
  else
  {
    // Found value.
    return TRUE;
  }
}

int
dynints_int_find(const int value, dynints_t *dynints)
{
  int i;

  for (i = 0; i < dynints->count; i++)
  {
    if (value == dynints->ints.a[i])
    {
      return i;
    }
  }

  // Didn't find value.
  return -1;
}

/**
** Check that the list of integers are the same.
** Each integer must be checked seperately becuase they may be in
** a different order.
**
** @param a_states first object to compare
** @param b_states second object to compare
** @return 0 when different, 1 when equivalent.
*/
boolean
dynints_compare(dynints_t *a_states, dynints_t *b_states)
{
  int a_states_count;
  int i;
  int *pcurrent_a_state;

  a_states_count = a_states->count;

  if (a_states_count != b_states->count)
  {
    // Don't have the same number of states, so they don't match.
    return FALSE;
  }

  for (i = 0, pcurrent_a_state = a_states->ints.a;
       i < a_states_count;
       i++, pcurrent_a_state++)
  {
    boolean status;

    // Search for state because they could be in a different order.
    status = dynints_int_exists(*pcurrent_a_state, b_states);
    if (FALSE == status)
    {
      // State not found
      return FALSE;
    }
  }

  return TRUE;
}

/**
** b_states should just contain one state.
** If a_states doesn't have anything, allow anything in b_states.
** If a_states just contains COLOR_STATE_WILDCARD, is true when any state
** is present in b_states.
** If a_states lists the state contained in b_states, is true.
** Otherwise, is false.
**
** @param a_states first object to compare
** @param b_states second object to compare
** @return 0 when different, 1 when equivalent.
*/
boolean
dynints_compare_cond(dynints_t *a_states, dynints_t *b_states)
{
  int b_state;
  int a_states_count;
  int i;
  int *pcurrent_a_state;

  if (0 == a_states->count)
  {
    return TRUE;
  }

  if (b_states->count > 1)
  {
    // Should never contain more than one state.
    error_printf("Found more than one state in a specie.\n");
    assert(0);
  }

  b_state = b_states->ints.a[0];

  a_states_count = a_states->count;

  for (i = 0, pcurrent_a_state = a_states->ints.a;
       i < a_states_count;
       i++, pcurrent_a_state++)
  {
    boolean status;

    // Search for state because they could be in a different order.
    status = dynints_int_exists(*pcurrent_a_state, b_states);
    if (*pcurrent_a_state == b_state)
    {
      // State found
      return TRUE;
    }
  }

  // b_state not listed in a_states.
  return FALSE;
}

boolean
dynints_int_remove(const int value, dynints_t *dynints)
{
  int i;
  int j;
  size_t count;
  int *ints;

  count = dynints->count;
  ints = dynints->ints.a;

  for (i = 0; i < count; i++)
  {
    if (value == ints[i])
    {
      for (j = i + 1; j < count; j++)
      {
        // Remove value by moving all values after current.
        ints[j - 1] = ints[j];
      }

      // Delete last element.
      ints[count - 1] = 0;

      // One element was removed.
      dynints->count--;

      return TRUE;
    }
  }

  // Didn't find value.
  return FALSE;
}

int *
dynints_int_alloc(dynints_t *dynints)
{
  int status;

  status = dynints_resize(dynints, dynints->count + 1);
  if (0 == status)
  {
    assert(0);
    return NULL;
  }

  return &dynints->ints.a[dynints->count++];
}

dynints_t *
dynints_create(const size_t init_count)
{
  dynints_t *new_dynints;
  int status;

  new_dynints = calloc(1, sizeof(dynints_t));
  if (NULL == new_dynints)
  {
    error_printf("Unable to allocate new dynints.\n");
    return NULL;
  }

  status = dynarray_create(init_count,
                           &new_dynints->alloc_count, sizeof(int),
                           &new_dynints->ints.v);
  if (0 == status)
  {
    // Error message output in dynarray_create.
    free(new_dynints);
    return NULL;
  }

  new_dynints->count = 0;

  return new_dynints;
}

boolean
dynints_resize(dynints_t *dynints, const size_t count)
{
  int status;

  status = dynarray_resize(count, &dynints->alloc_count, sizeof(int),
                           &dynints->ints.v);
  if (0 == status)
  {
    // Error message output in dynarray_resize.
    return FALSE;
  }

  return TRUE;
}

void
dynints_destroy(dynints_t *dynints)
{
  free(dynints->ints.a);
  free(dynints);
}

