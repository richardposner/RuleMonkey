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

#ifndef DYNRULES_H
#define DYNRULES_H

typedef struct dynrules_struct dynrules_t;

#define DYNRULES_OUTPUT(ptype, level, ...) \
   CHECKDEBUG(dynrules_output, ptype, level, __VA_ARGS__)

#define DYNRULES_OUTPUT_RULE(ptype, level, ...) \
   CHECKDEBUG(dynrules_output_rule, ptype, level, __VA_ARGS__)

typedef union rulearray_union rulearray_t;

#include "rule.h"
#include "dynmaps.h"
#include "dynconnections.h"
#include "dyncolors.h"
#include "dynints.h"
#include "output.h"

union rulearray_union
{
  rule_t *a;
  void *v;
};

struct dynrules_struct
{
  // Pointer to allocated array of rules.
  rulearray_t rules;

  // Number of elements currently in use.
  size_t count;

  // Number of elements currently allocated.
  size_t alloc_count;
};

int dynrules_same_components
      (rule_t *rule_a, rule_t *rule_b, dyncolors_t *dyncolors);
boolean dynrules_rule_create(const int color_value,
                             int modifier_flags,
                             const count_flag_t count_flag, const int count,
                             const double delta,
                             const int xor_index,
                             const int or_index,
                             const int adjacent_index,
                             const int and_index,
                             dynrules_t *dynrules);
boolean dynrules_rule_add_modifier_flags(const int index, const int or_mask,
                                         dynrules_t *dynrules);
boolean dynrules_rule_remove_modifier_flags(const int index,
                                            const int and_mask,
                                            dynrules_t *dynrules);
boolean dynrules_rule_modify_count(const int index,
                                   const count_flag_t count_flag,
                                   const int count, dynrules_t *dynrules);
boolean dynrules_rule_connect_any_one_rule(const int index,
                                           dynrules_t *dynrules);
boolean dynrules_rule_connect_zero_or_one_rule(const int index,
                                               dynrules_t *dynrules);
void dynrules_rule_add_xor_rule(const int index, const int xor_index,
                                dynrules_t *dynrules);
void dynrules_rule_add_or_rule(const int index, const int or_index,
                               dynrules_t *dynrules);
void dynrules_rule_add_adjacent_rule(const int index,
                                     const int adjacent_index,
                                     dynrules_t *dynrules);
void dynrules_rule_add_and_rule(const int index, const int and_index,
                                 dynrules_t *dynrules);

boolean dynrules_rule_remove_xor_rule(const int index, const int and_index,
                                      dynrules_t *dynrules);
boolean dynrules_rule_remove_or_rule(const int index, const int and_index,
                                     dynrules_t *dynrules);
boolean dynrules_rule_remove_adjacent_rule(const int index,
                                           const int adjacent_index,
                                           dynrules_t *dynrules);
boolean dynrules_rule_remove_and_rule(const int index, const int and_index,
                                       dynrules_t *dynrules);

rule_t *dynrules_rule_alloc(dynrules_t *dynrules);
int dynrules_rule_add_bind_index(const int index, const int bind_index,
                                 dynrules_t *dynrules);
boolean dynrules_rule_add_primary_self_bind_index(const int index,
                                                  const int bind_index,
                                                  dynrules_t *dynrules);
boolean dynrules_rule_add_secondary_self_bind_index(const int index,
                                                    const int bind_index,
                                                    dynrules_t *dynrules);
boolean dynrules_rule_add_primary_split_index(const int index,
                                              const int bind_index,
                                              dynrules_t *dynrules);
boolean dynrules_rule_add_secondary_split_index(const int index,
                                                const int bind_index,
                                                dynrules_t *dynrules);
boolean dynrules_map_molecules(const int lhs_count, const int rhs_count,
                               const int molecule_star_lhs,
                               const int molecule_star_rhs,
                               int *specie_rule_indexes,
                               dynmaps_t *dynmaps,
                               dynints_t *extra_lhs_molecules,
                               dynints_t *extra_rhs_molecules,
                               dynconnections_t *dynconnections,
                               dyncolors_t *dyncolors, dynrules_t *dynrules);
boolean dynrules_find_matching_molecule
          (size_t current_lhs_rule_index,
           size_t rhs_index, dynmaps_t *dynmaps,
           dynconnections_t *dynconnections,
           dyncolors_t *dyncolors, dynrules_t *dynrules);
boolean dynrules_compare_rules
          (size_t lhs_rule_index, size_t rhs_rule_index,
           size_t rhs_index,
           dynconnections_t *dynconnections,
           dyncolors_t *dyncolors, dynrules_t *dynrules);
boolean dynrules_complete(dynrules_t *dynrules);
dynrules_t *dynrules_create(const size_t init_count);
void dynrules_output_rule(const output_type ptype, const int level,
                          rule_t *rule, dynrules_t *dynrules,
                          dyncolors_t *dyncolors);
void dynrules_output(const output_type ptype, const int level,
                     dynrules_t *dynrules, dyncolors_t *dyncolors);
void dynrules_destroy(dynrules_t *da);

#endif
