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


#include "fsa.h"
#include "output.h"
#include "constants.h"
#include "dynarray.h"
#include "llist.h"
#include <assert.h>

void suppress_warnings1() {SUPPRESS_WARNINGS;}

fsa_t *
fsa_create(size_t elem_size, size_t elem_per_block,
           fsa_callback_t init_callback, fsa_callback_t clean_callback)
{
  fsa_t *fsa;
  int status;

  fsa = calloc(1, sizeof(fsa_t));

  fsa->elem_size = elem_size;
  fsa->elem_per_block = elem_per_block;
  fsa->last_count = elem_per_block; // So that new block will be created.

  fsa->init_callback = init_callback;
  fsa->clean_callback = clean_callback;

  status = dynarray_create(INITIAL_FSA_BLOCK_COUNT, &fsa->block_alloc_count,
                           sizeof(void *), &fsa->blocks.v);
  if (FALSE == status)
  {
    assert(0);
    free(fsa);
    return NULL;
  }

  return fsa;
}

void
fsa_destroy(fsa_t *fsa)
{
  int i;

  if (fsa->clean_callback != NULL)
  {
    fsa_traverse(fsa->clean_callback, NULL, fsa);
  }

  llist_destroy(fsa->available);
  for (i = 0; i < fsa->block_count; i++)
  {
    free(fsa->blocks.a[i]);
  }
  free(fsa->blocks.v);
  free(fsa);
}

void *
fsa_get_new(fsa_t *fsa)
{
  void *pelem;
  int index;

  if (fsa->available != NULL)
  {
    // Previously removed element is available for reuse.
    pelem = llist_shift(&fsa->available);

    // Don't know index of reused element.  If important, should be
    // preserved inside data element.
    index = -1;
  }
  else if(fsa->elem_per_block == fsa->last_count)
  {
    int status;

    // No elements available, need to create new block.
    fsa->block_count++;
    status = dynarray_resize(fsa->block_count, &fsa->block_alloc_count,
                             sizeof(void *), &fsa->blocks.v);
    if (FALSE == status) return NULL;
    pelem = calloc(fsa->elem_size, fsa->elem_per_block);
    if (NULL == pelem) return NULL;
    fsa->blocks.a[fsa->block_count - 1] = pelem;
    fsa->plast_elem = pelem;
    fsa->last_count = 1;

    index =
      ((fsa->block_count - 1) * fsa->elem_per_block) + fsa->last_count - 1;
  }
  else
  {
    // Take next available element in last block.
    pelem = fsa->plast_elem + fsa->elem_size;
    fsa->plast_elem = pelem;
    fsa->last_count++;

    index =
      ((fsa->block_count - 1) * fsa->elem_per_block) + fsa->last_count - 1;
  }

  if (fsa->init_callback != NULL)
  {
    fsa->init_callback(&index, pelem);
  }

  // Return pointer to new element.
  return pelem;
}

int
fsa_remove(void *elem, fsa_t *fsa)
{
  int status;

  if (fsa->clean_callback != NULL)
  {
    fsa->clean_callback(NULL, elem);
  }

  status = llist_push(elem, &fsa->available);

  return status;
}

size_t
fsa_count(fsa_t *fsa)
{
  size_t count;

  count =  (fsa->elem_per_block * fsa->block_count) -
           (fsa->elem_per_block - fsa->last_count) -
           llist_count(fsa->available);

  return count;
}

void *
fsa_get_existing(size_t index, fsa_t *fsa)
{
  size_t block_index;
  size_t elem_index;
  char *block; // Use (char *) so pointer arithmetic works properly.

  block_index = index / fsa->elem_per_block;
  assert(block_index <= fsa->block_count);
  elem_index = index - (fsa->elem_per_block * block_index);

  block = fsa->blocks.a[block_index];
  
  return block + elem_index * fsa->elem_size;
}

/**
 * Calls fsa_callback with all elements.
 * Note this also includes removed elements.
 */
int
fsa_traverse(fsa_callback_t fsa_callback, void *d, fsa_t *fsa)
{
  size_t i_block;
  char **pcurrent_block;
  size_t i_elem;
  void *current_elem;
  int elem_count;

  for (i_block = 0, pcurrent_block = fsa->blocks.a;
       i_block < fsa->block_count;
       i_block++, pcurrent_block++)
  {
    if (i_block == (fsa->block_count - 1))
    {
      // Last block, only go through currently used elements.
      elem_count = fsa->last_count;
    }
    else
    {
      // Not the last block, go through all used elements.
      elem_count = fsa->elem_per_block;
    }

    for (i_elem = 0, current_elem = *pcurrent_block;
         i_elem < elem_count;
         i_elem += 1, current_elem += fsa->elem_size)
    {
      int status;

      status = fsa_callback(d, current_elem);
      if (1 != status)
      {
        return status;
      }
    }
  }

  return 1;
}
