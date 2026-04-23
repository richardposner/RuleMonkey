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

#ifndef FSA_H
#define FSA_H

#include "llist.h"
#include <stdio.h>

#ifndef NDEBUG
#define FSA_CREATE(es, elem_per_block, init_callback, clean_callback) \
  fsa_create((es), elem_per_block, init_callback, clean_callback)
#define FSA_GET_NEW(es, fsa) \
  ((es) != (fsa)->elem_size) ? assert(0),abort(),NULL \
  : fsa_get_new(fsa) 
#define FSA_DESTROY(es, fsa)  \
  ((es) != (fsa)->elem_size) ? assert(0),abort(),1 \
  : fsa_destroy(fsa) 
#define FSA_REMOVE(es, elem, fsa)  \
  ((es) != (fsa)->elem_size) ? assert(0),abort(),1 \
  : fsa_remove(elem, fsa) 
#define FSA_COUNT(es, fsa)  \
  ((es) != (fsa)->elem_size) ? assert(0),abort(),1 \
  : fsa_count(fsa) 
#define FSA_TRAVERSE(es, fsa_callback, d, fsa)  \
  ((es) != (fsa)->elem_size) ? assert(0),abort(),1 \
  : fsa_traverse(fsa_callback, d, fsa) 
#define FSA_GET_EXISTING(es, index, fsa)  \
  ((es) != (fsa)->elem_size) ? assert(0),abort(),NULL \
  : fsa_get_existing(index, fsa) 
#else
#define FSA_CREATE(es, elem_per_block, init_callback, clean_callback) \
  fsa_create(es, elem_per_block, init_callback, clean_callback)
#define FSA_GET_NEW(es, fsa) fsa_get_new(fsa)
#define FSA_DESTROY(es, fsa) fsa_destroy(fsa)
#define FSA_REMOVE(es, elem, fsa) fsa_remove(elem, fsa)
#define FSA_COUNT(es, fsa) fsa_count(fsa)
#define FSA_TRAVERSE(es, fsa_callback, d, fsa) \
  fsa_traverse(fsa_callback, d, fsa)
#define FSA_GET_EXISTING(es, index, fsa) fsa_get_existing(index, fsa)
#endif

typedef struct fsa_struct fsa_t;
typedef union fsa_blocks_union fsa_blocks_t;

typedef int (*fsa_callback_t) (void *d, void *elem);

union fsa_blocks_union
{
  char **a;
  void *v;
};


struct fsa_struct
{
  // Size of each element.
  size_t elem_size;

  // Number of elements in each block.
  size_t elem_per_block;

  // Number of blocks currently in use.
  size_t block_count;

  // Number of blocks currently allocated.
  size_t block_alloc_count;

  // Number of elements used in last block.
  size_t last_count;

  // Pointer to last new element returned from block.
  // This does not get updated when reused element is returned.
  void *plast_elem;

  // Function to initialize any new memory allocations before returning
  // the pointer.
  // Note that the first parameter will always get set to NULL.
  fsa_callback_t init_callback;

  // Function to cleanup any memory allocations before removing element
  // or destroying complete fsa structure.
  // Note that the first parameter will always get set to NULL.
  fsa_callback_t clean_callback;

  // List of allocated but currently unused entries.
  llist_t *available;

  // Dynamic array of data blocks.
  // Size of each block = (elem_size * elem_per_block).
  // Note that array may be reallocated, but blocks are not.
  // This means that it is safe to store a pointer to any given element.
  // Use char ** so that pointer arithmetic is handled properly.
  fsa_blocks_t blocks;
};


fsa_t *fsa_create(size_t elem_size, size_t elem_per_block,
                  fsa_callback_t init_callback, fsa_callback_t clean_callback);
void fsa_destroy(fsa_t *fsa);
void *fsa_get_new(fsa_t *fsa);
int fsa_remove(void *elem, fsa_t *fsa);
size_t fsa_count(fsa_t *fsa);
int fsa_traverse(fsa_callback_t fsa_callback, void *d, fsa_t *fsa);
void *fsa_get_existing(size_t index, fsa_t *fsa);

#endif
