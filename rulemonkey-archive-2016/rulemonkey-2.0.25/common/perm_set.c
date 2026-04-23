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


#include "perm_set.h"
#include "fsa.h"
#include "constants.h"
#include "dynints.h"
#include <assert.h>

void suppress_warnings2() {SUPPRESS_WARNINGS;}

int
perm_set_unique_check(const int set_size, int *current_perm)
{ 
  int i;
  int current_node;
  
  for (i = 0; i < set_size; i++)
  {
    int j;
  
    current_node = current_perm[i];
  
    for (j = i + 1; j < set_size; j++)
    {
      if (current_node == current_perm[j])
      {
        // Single node assigned to more than on component.
        return FALSE;
      }                
    }
  }
  
  // Each node is only assigned to one component.
  return TRUE;
}   

/**
 * Create new perm_set to collect all component node permutations.
 */
perm_set_t *
perm_set_create(const int set_size)
{
  perm_set_t *new_perm_set;

  new_perm_set = calloc(sizeof(perm_set_t), 1);
  if (NULL == new_perm_set) return NULL;

  new_perm_set->set_size = set_size;
  new_perm_set->elem_size = set_size * sizeof(int);

  new_perm_set->fsa = fsa_create(new_perm_set->elem_size,
                                  INITIAL_COMP_PERM_COUNT, NULL, NULL);
  if (NULL == new_perm_set->fsa)
  {
    free(new_perm_set);
    return NULL;
  }

  return new_perm_set;
}

void
perm_set_output(const output_type ptype, const int level, perm_set_t *perm_set)
{
  int i;
  int j;
  int count;

  count = FSA_COUNT(perm_set->elem_size, perm_set->fsa);

  for (i = 0; i < count; i++)
  {
    int *current_perm;

    current_perm = fsa_get_existing(i, perm_set->fsa);

    CUSTOM_PRINTF(ptype, level, "[");
    for (j = 0; j < perm_set->set_size; j++)
    {
      if (j != 0)
      {
        CUSTOM_PRINTF(ptype, level, ",");
      }
      CUSTOM_PRINTF(ptype, level, "%d", current_perm[j]);
    }
    CUSTOM_PRINTF(ptype, level, "]");
  }
}

void
perm_set_destroy(perm_set_t *perm_set)
{
  FSA_DESTROY(perm_set->elem_size, perm_set->fsa);
  free(perm_set);
}

/**
 * Permute all node_map so that every combination of nodes is produced.
 * Note that this may create impossible states (single node assigned to more
 * than one component).  If we were to fix this now, it would be very
 * difficult to ensure that all permutations were generated.  This means that
 * impossible states will be removed later.
 */
int
perm_set_add_dynints(const int perm_base, const int molecule_node,
                     const int index, dynints_t *node_map,
                     perm_set_t *perm_set)
{
  int i;
  int j;
  int *pnode;
  int first_node;
  int count;

  // Note that extra permutations will be added while executing this function,
  // but we don't need to process new entries, so count is not updated.
  count = FSA_COUNT(perm_set->elem_size, perm_set->fsa) - perm_base;

  assert (count >= 0);

  if (count == 0)
  {
    int *current_perm;

    // Create first permutation.
    current_perm = FSA_GET_NEW(perm_set->elem_size, perm_set->fsa);
    count = FSA_COUNT(perm_set->elem_size, perm_set->fsa) - perm_base;
    if (count != 1) return FALSE;

    // Set molecule node in first permutation.
    current_perm[0] = molecule_node;
  }

  if (NULL == node_map)
  {
    // No components present, so nothing more to do.
    return TRUE;
  }

  assert(node_map->count > 0);

  pnode = node_map->ints.a;
  first_node = *pnode;

  // First node doesn't need to be permuted, so simply add to all permutations.
  for (i = 0; i < count; i++)
  {
    int *current_perm;

    current_perm = FSA_GET_EXISTING(perm_set->elem_size, i + perm_base,
                                    perm_set->fsa);
    current_perm[index] = first_node;
    assert(molecule_node == current_perm[0]);
  }

  // Loop through each node mapping, we've already taken care of the first one.
  for (j = 1, pnode++;
       j < node_map->count;
       j++, pnode++)
  {
    int k;

    // Loop through each permutation.  Note that new permutations will be
    // created while this loop is executing, but we only need to go through
    // the original permutations, so count is not upated.
    for (k = 0; k < count; k++)
    {
      int *current_perm;

      current_perm = FSA_GET_EXISTING(perm_set->elem_size, k + perm_base,
                                      perm_set->fsa);

      // Permute off of first node.
      if (current_perm[index] == first_node)
      {
        int *new_perm;

        new_perm = FSA_GET_NEW(perm_set->elem_size, perm_set->fsa);
        if (NULL == new_perm)
        {
          assert(0);
          return FALSE;
        }

        // Copy existing permutation to new permutation.
        memcpy(new_perm, current_perm, perm_set->elem_size);

        // Replace with new node.
        new_perm[index] = *pnode;
      }
    }
  }

  return TRUE;
}

int
perm_set_count(perm_set_t *perm_set)
{
  return FSA_COUNT(perm_set->elem_size, perm_set->fsa);
}

int *
perm_set_get_existing(const int index, perm_set_t *perm_set)
{
  return FSA_GET_EXISTING(perm_set->elem_size, index, perm_set->fsa);
}

