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

/*
 * Handles generating all component rule -> node map permutations.
 */
#ifndef PERM_SET_H
#define PERM_SET_H

#define PERM_SET_OUTPUT(ptype, level, ...) \
  CHECKDEBUG(perm_set_output, ptype, level, __VA_ARGS__)

typedef struct perm_set_struct perm_set_t;

#include "rule.h"
#include "fsa.h"
#include "output.h"
#include "dynints.h"

struct perm_set_struct
{
  // Width of array element (components in current molecule + molecule).
  int set_size;

  // Size of array element (comp_count * sizeof(int)).
  size_t elem_size;

  // All permutations.
  fsa_t *fsa;

  // Pointer to currently select permutation.
  int *current_perm;
};

int perm_set_generate_all(rule_t *molecule_rule, const int molecule_node);
void perm_set_destroy(perm_set_t *perm_set);
void perm_set_output(const output_type ptype, const int level,
                      perm_set_t *perm_set);
int *perm_set_get_existing(const int index, perm_set_t *perm_set);
int perm_set_count(perm_set_t *perm_set);
int perm_set_unique_check(const int set_size, int *current_perm);
perm_set_t *perm_set_create(const int set_size);
int perm_set_add_dynints(const int perm_base, const int molecule_node,
                         const int index, dynints_t *node_map,
                         perm_set_t *perm_set);

#endif
