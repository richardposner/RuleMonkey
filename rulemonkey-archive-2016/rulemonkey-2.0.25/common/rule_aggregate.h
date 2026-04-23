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

#ifndef RULE_AGGREGATE_H
#define RULE_AGGREGATE_H

#include "cgraph.h"
#include "rule.h"
#include "output.h"

#define AGG_MULTIPLIER_OFFSET 0
#define AGG_USED_OFFSET 0
#define AGG_NODE_OFFSET 1
#define AGG_PRI_USED_OFFSET 0
#define AGG_PRI_NODE_OFFSET 1
#define AGG_SEC_USED_OFFSET 2
#define AGG_SEC_NODE_OFFSET 3
#define AGG_OLD_STATE_OFFSET 2
#define AGG_NEW_STATE_OFFSET 3

#define AGG_SINGLE_LENGTH 2
#define AGG_DOUBLE_LENGTH 4
#define FULL_AGG_SINGLE_LENGTH 2
#define FULL_AGG_DOUBLE_LENGTH 4
#define AGG_BIND_LENGTH AGG_SINGLE_LENGTH
#define FULL_AGG_BIND_LENGTH FULL_AGG_DOUBLE_LENGTH
#define AGG_SPLIT_LENGTH AGG_DOUBLE_LENGTH
#define FULL_AGG_SPLIT_LENGTH FULL_AGG_DOUBLE_LENGTH
#define AGG_SELF_BIND_LENGTH AGG_DOUBLE_LENGTH
#define FULL_AGG_SELF_BIND_LENGTH FULL_AGG_DOUBLE_LENGTH
#define AGG_STATE_CHANGE_LENGTH 4
#define FULL_AGG_STATE_CHANGE_LENGTH 4

#define RULE_AGGREGATE_OUTPUT(ptype, level, full_aggregate) \
   CHECKDEBUG(rule_aggregate_output, ptype, level, full_aggregate)

#define RULE_AGGREGATE_OUTPUT_FULL(ptype, level, full_aggregate) \
   CHECKDEBUG(rule_aggregate_output_full, ptype, level, full_aggregate)


size_t
rule_aggregate_num_ints(const size_t split_count,
                        const size_t bind_count,
                        const size_t self_bind_count,
                        const size_t state_change_count,
                        const size_t transient_bind_count,
                        const size_t transient_self_bind_count);

void
rule_aggregate_init_full(int *full_rule_aggregate, int *rule_aggregate);

size_t
rule_aggregate_sizeof_full(const size_t split_count,
                           const size_t bind_count,
                           const size_t self_bind_count,
                           const size_t state_change_count,
                           const size_t transient_bind_count,
                           const size_t transient_self_bind_count);

boolean
rule_aggregate_cgraph_append(cgraph_t *a, int *rule_aggregate_full,
                             cgraph_t *b, int *rule_aggregate_b,
                             cgraph_t *c);

boolean
rule_aggregate_cgraph_append_raw(cgraph_t *a, int *rule_aggregate_full,
                                 cgraph_t *b, int *rule_aggregate_b,
                                 cgraph_t *c, int *renum_a, int *renum_b);

boolean
rule_aggregate_cgraph_apply(cgraph_t *a, int *rule_aggregate_full);

boolean
rule_aggregate_compare(const int *aggregate1, const int *aggregate2);

boolean
rule_aggregate_first_bind(int *aggregate, int *ppri);

boolean
rule_aggregate_first_active(int *aggregate, int *ppri, int *psec, int *pid,
                            boolean *pbimolecular);

void
rule_aggregate_output_raw(const output_type ptype, const int level,
                          int *rule_aggregate, size_t num_ints);

void
rule_aggregate_output(const output_type ptype, const int level,
                      int *rule_aggregate);

void
rule_aggregate_output_full(const output_type ptype, const int level,
                           int *rule_aggregate);

int *
rule_aggregate_state_change_add_rule(int *rule_aggregate,
                                     const size_t count,
                                     ilist_t *current_node,
                                     ilist_t *old_list,
                                     ilist_t *new_list,
                                     const int matched_node, int *pconflict);

int *
rule_aggregate_single_add_rule(int *rule_aggregate,
                               const size_t count,
                               ilist_t *current_node,
                               const int matched_node, int *changed);

int *
rule_aggregate_double_add_rule(int *rule_aggregate,
                               const size_t count,
                               const ilist_t *current_pri_node,
                               const ilist_t *current_sec_node,
                               const int matched_node, int *changed);

int rule_aggregate_calculate_product_count(cgraph_t *a, int *rule_aggregate);

#endif
