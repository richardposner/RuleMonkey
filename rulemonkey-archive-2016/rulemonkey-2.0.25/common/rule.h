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

#ifndef RULE_H
#define RULE_H

typedef enum count_flag_enum count_flag_t;
typedef struct rule_struct rule_t;

#include "perm_set.h"
#include "specie.h"
#include "dynrules.h"
#include "../nauty22/nauty.h"
#include "ilist.h"
#include "rule.h"
#include "dynrules.h"
#include "dynints.h"
#include "dynrule_aggregates.h"
#include "constants.h"
#include "output.h"
#include "dyncolors.h"

// Special value used to match any color.
#define SPECIE_Any -1

// Define bitmasks to modify behavior of rules and/or reactants.
#define NO_MODIFIER_FLAGS           0x0000
#define ADJACENT_EXCLUSIVE_MODIFIER 0x0001

#define MAX_COUNT_FLAG_STR_LENGTH (sizeof("GREAT_EQ") + 1)
enum count_flag_enum
{
  IGNORE = 0,
  EQUAL,
  LESS_EQ,
  GREAT_EQ,
};

#define RULE_OUTPUT(ptype, level, ...) \
  CHECKDEBUG(rule_output, ptype, level, __VA_ARGS__)


struct rule_struct
{
  /* Attributes that node must meet to match rule */

  // Index of rule.
  int index;

  // Color rule is looking for.
  int color_value;

  // List of states rule is looking for.
  dynints_t *state_dynints;

  // Assorted flags to modify behavior of rules and/or reactants.
  int modifier_flags;

  // Modifier for count.
  count_flag_t count_flag;

  // Number of matching nodes that must be present.
  // Undefined if count_flag == IGNORE.
  int count;

  // Adjacent rule.
  int adjacent_index;

  // Pointer to adjacent rule.
  rule_t *adjacent_rule;

  // And rule.
  int and_index;

  // Pointer to and rule.
  rule_t *and_rule;

  // Or rule.
  int or_index;

  // Pointer to or rule.
  rule_t *or_rule;

  // Exclusive or rule.
  int xor_index;

  // Pointer to exclusive or rule.
  rule_t *xor_rule;

  /* Attributes to apply to matching node */

  // TRUE if rule contains attributes that modify node.
  int attributes_modified;

  // List of indexes to match up other nodes to split with.
  ilist_t *primary_split_index_list;
  ilist_t *secondary_split_index_list;

  // List of indexes to match up other nodes to bind with.
  ilist_t *bind_index_list;

  // List of indexes to match up other nodes to split with.
  ilist_t *primary_self_bind_index_list;
  ilist_t *secondary_self_bind_index_list;

  // List of indexes to ensure accurate permutations for state changes.
  // Map of existing states and what they need to be changed to.
  // Note that all entries in state_change_old_list are expected to be
  // listed in state_dynints above.
  ilist_t *state_change_index_list;
  ilist_t *state_change_old_list;
  ilist_t *state_change_new_list;

  // List of indexes to match up other nodes to bind with.
  ilist_t *transient_bind_index_list;

  // List of indexes to match up other nodes to split with.
  ilist_t *primary_transient_self_bind_index_list;
  ilist_t *secondary_transient_self_bind_index_list;

  // Used to help match nodes to rules.
  dynints_t *node_map;

  // Component permutations associated with current molecule.
  perm_set_t *perm_set;

  // Index into perm_set for current component or molecule
  // If molecule rule has no components, contains molecule node index.
  int perm_set_index;
};

int rule_check(specie_t *specie, dynrule_aggregates_t *dynrule_aggregates,
               rule_t *rule, int *pcount,
               dynrules_t *dynrules, dyncolors_t *dyncolors);
boolean rule_setup(const int index, const int color_value, int modifier_flags,
                   const count_flag_t count_flag, const int count,
                   const double delta,
                   const int xor_index, const int or_index,
                   const int adjacent_index, const int and_index,
                   rule_t *new_rule);
void rule_output(const output_type ptype, const int level,
                 rule_t *rule, dyncolors_t *dyncolors,
                 const char *prefix, boolean *used);
int rule_add_bind_index(const int index, rule_t *rule);
int rule_add_primary_self_bind_index(const int index, rule_t *rule);
int rule_add_secondary_self_bind_index(const int index, rule_t *rule);
int rule_add_state_change_index(const int index, rule_t *rule);
int rule_add_primary_split_index(const int index, rule_t *rule);
int rule_add_secondary_split_index(const int index, rule_t *rule);
int rule_add_state(const int state, const int new_state, rule_t *rule);
int rule_add_modify_state(const int index, const int state,
                          const int new_state, rule_t *rule);
int rule_get_max_primary_split_elements(rule_t *rule, char *used_rules);
int rule_get_max_secondary_split_elements(rule_t *rule, char *used_rules);
int rule_get_max_bind_elements(rule_t *rule, char *used_rules);
int rule_get_max_primary_self_bind_elements(rule_t *rule, char *used_rules);
int rule_get_max_secondary_self_bind_elements(rule_t *rule, char *used_rules);
int rule_get_max_state_change_elements(rule_t *rule, char *used_rules);
int rule_get_max_transient_bind_elements(rule_t *rule, char *used_rules);
int rule_get_max_primary_transient_self_bind_elements(rule_t *rule,
                                                      char *used_rules);
int rule_get_max_secondary_transient_self_bind_elements(rule_t *rule,
                                                        char *used_rules);
int rule_create_specie(rule_t *rule, int *pnew_specie_id,
                       dynspecies_t *dynspecies,
                       dynrules_t *dynrules);
int rule_complete(rule_t *rule, dynrules_t *dynrules);
void rule_destroy(rule_t *rule);


#endif
