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

#ifndef DYNRULE_AGGREGATES_H
#define DYNRULE_AGGREGATES_H

typedef struct basic_aggregate_struct basic_aggregate_t;
typedef struct dynrule_aggregates_struct dynrule_aggregates_t;

#include "dynrules.h"
#include "rule.h"
#include "specie.h"
#include "reaction_definition.h"
#include "output.h"

#define DYNRULE_AGGREGATES_OUTPUT(ptype, level, ...) \
   CHECKDEBUG(dynrule_aggregates_output, ptype, level, __VA_ARGS__)

struct dynrule_aggregates_struct
{
  size_t split_count;
  size_t bind_count;
  size_t self_bind_count;
  size_t state_change_count;
  size_t transient_bind_count;
  size_t transient_self_bind_count;

  /*
   * (1 +
   *  (split_count * AGG_SPLIT_LENGTH) + 1 +
   *  (bind_count * AGG_BIND_LENGTH) + 1 +
   *  (self_bind_count * AGG_SELF_BIND_LENGTH) + 1 +
   *  (state_change_count * STATE_CHANGE_COUNT_LENGTH) + 1 +
   *  (transient_bind_count * AGG_BIND_LENGTH) + 1 +
   *  (transient_self_bind_count * AGG_SELF_BIND_LENGTH) + 1
   * ) * sizeof(int)
   *
   *
   * multiplier
   *
   * split block:
   *   int primary_used;        // Offset 0
   *   int primary_matched_node;// Offset 1
   *   int secondary_used;      // Offset 2
   *   int secondary_matched_node; // Offset 3
   *   ...
   * int terminator = -1;
   *
   * bind block:
   *   int used;                // Offset 0
   *   int matched_node;        // Offset 1
   *   ...
   * int terminator = -1;
   *
   * self_bind block:
   *   int primary_used;        // Offset 0
   *   int primary_matched_node;// Offset 1
   *   int secondary_used;      // Offset 2
   *   int secondary_matched_node; // Offset 3
   *   ...
   * int terminator = -1;
   *
   * state change block:
   *   int used;                // Offset 0
   *   int matched_node;        // Offset 1
   *   int old_color;           // Offset 2
   *   int new_color;           // Offset 3
   *   ...
   * int terminator = -1;
   *
   * transient bind block:
   *   int used;                // Offset 0
   *   int matched_node;        // Offset 1
   *   ...
   * int terminator = -1;
   *
   * transient_self_bind block:
   *   int primary_used;        // Offset 0
   *   int primary_matched_node;// Offset 1
   *   int secondary_used;      // Offset 2
   *   int secondary_matched_node; // Offset 3
   *   ...
   * int terminator = -1;
   *
   */
  size_t rule_aggregate_num_ints;
  size_t rule_aggregate_size;

  size_t count;
  size_t alloc_count;

  // All rule permutations must produce the correct number of products.
  int product_count;

  // contains a single aggregate that is used to initialize all new aggregates.
  int *rule_aggregate_template;

  // Block of all aggregates.
  intarray_t rule_aggregates;

  // Pointer to last aggregate to be allocated.
  int *current_aggregate;
};

int dynrule_aggregates_identical
      (dynrule_aggregates_t *a, dynrule_aggregates_t *b);

dynrule_aggregates_t *dynrule_aggregates_create
                    (const size_t init_count,
                     const size_t split_count,
                     const size_t bind_count,
                     const size_t self_bind_count,
                     const size_t state_change_count,
                     const size_t transient_bind_count,
                     const size_t transient_self_bind_count,
                     const int product_count);
void dynrule_aggregates_destroy(dynrule_aggregates_t *dynrule_aggregates);
size_t dynrule_aggregates_sizeof_full_rule_aggregate
         (dynrule_aggregates_t *dynrule_aggregates);
void dynrule_aggregates_reset_all(dynrule_aggregates_t *dynrule_aggregates);
void dynrule_aggregates_output(const output_type ptype, const int level,
                               dynrule_aggregates_t *dynrule_aggregates);
int dynrule_aggregates_add_rule (dynrule_aggregates_t *dynrule_aggregates,
                                 rule_t *rule, const int matched_node);
int dynrule_aggregates_check_current_unique
      (dynrule_aggregates_t *dynrule_aggregates);

int dynrule_aggregates_check_product_count
      (dynrule_aggregates_t *dynrule_aggregates, specie_t *specie);

int dynrule_aggregates_rule_aggregate_create
        (dynrule_aggregates_t *dynrule_aggregates);
int dynrule_aggregates_count_shared_components
      (dynrule_aggregates_t *a, dynrule_aggregates_t *b,
       specie_t *specie, dyncolors_t *dyncolors,
       double *a_onlyp, double *b_onlyp, double *ab_bothp);

#endif
