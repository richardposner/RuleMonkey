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


#include "reaction_definition.h"
#include <assert.h>
#include "constants.h"
#include "output.h"

void suppress_warnings8() {SUPPRESS_WARNINGS;}

static int
reaction_definition_complete_raw(reaction_definition_t *reaction_definition,
                                 dynrules_t *dynrules,
                                 const size_t used_rules_size,
                                 char *used_rules);

static reaction_definition_t *
reaction_definition_create(const char *name, const double rate_const,
                           const int flags, const int reactant_count,
                           const int product_count,
                           dynints_t *new_species,
                           reaction_definition_t *next);

int
reaction_definition_add(const char *name, const double rate_const,
                        const int flags, const int reactant_count,
                        const int product_count,
                        dynints_t *new_species,
                        reaction_definition_t **pnext)
{
  reaction_definition_t *new_reaction_definition;

  new_reaction_definition =
    reaction_definition_create(name, rate_const, flags, reactant_count,
                               product_count, new_species, *pnext);

  if (NULL == new_reaction_definition)
  {
    assert(0);
    return 0;
  }

  *pnext = new_reaction_definition;

  return 1;
}

/*
** Returns pointer to newly allocated and initialized reaction definition.
*/
reaction_definition_t *
reaction_definition_create(const char *name, const double rate_const,
                           const int flags, const int reactant_count,
                           const int product_count,
                           dynints_t *new_species,
                           reaction_definition_t *next)
{
  reaction_definition_t *new_reaction_definition;

  new_reaction_definition = calloc((size_t)1, sizeof(reaction_definition_t));
  if (NULL == new_reaction_definition)
  {
    perror("Allocating new reaction definition");
    assert(0);
    return NULL;
  }

  // Set structure members.
  if (NULL == next)
  {
    new_reaction_definition->index = 1;
  }
  else
  {
    new_reaction_definition->index = next->index + 1;
  }
  new_reaction_definition->name = name ? strdup(name) : NULL;
  new_reaction_definition->rate_const = rate_const;
  new_reaction_definition->flags = flags;
  new_reaction_definition->reactant_definition_count = reactant_count;
  new_reaction_definition->product_count = product_count;
  new_reaction_definition->reactant_definitions =
    calloc((size_t)reactant_count, sizeof(reactant_definition_t));
  new_reaction_definition->reactant_llists =
    calloc((size_t)reactant_count, sizeof(llist_t *));
  new_reaction_definition->reactant_work_llists =
    calloc((size_t)reactant_count, sizeof(llist_t *));
  new_reaction_definition->reactant_work_flags =
    calloc((size_t)reactant_count, sizeof(int));
  new_reaction_definition->new_species = new_species;
  new_reaction_definition->next = next;

  return new_reaction_definition;
}

void
reaction_definition_output(const output_type ptype, const int level,
                           reaction_definition_t *reaction_definition,
                           dynrules_t *dynrules, dyncolors_t *dyncolors)
{
  int i;

  custom_printf(ptype, level,
                "Reaction definition: %s[%d](0x%.2x)\n",
                reaction_definition->name,
                reaction_definition->index, reaction_definition->flags);
  custom_printf(ptype, level, "Rate constant: %e\n",
                reaction_definition->rate_const);
  custom_printf(ptype, level, "Probability: %e\n", reaction_definition->prob);
  custom_printf
        (ptype, level,
         "Counts: s-%d, b-%d, sb-%d, sc-%d, tb-%d, tsb-%d\n",
         reaction_definition->split_count, reaction_definition->bind_count,
         reaction_definition->self_bind_count,
         reaction_definition->state_change_count,
         reaction_definition->transient_bind_count,
         reaction_definition->transient_self_bind_count);
  for (i = 0; i < reaction_definition->reactant_definition_count; i++)
  {
    custom_printf(ptype, level, "Reactant %d:", i + 1);
    custom_printf(ptype, level, "\n");

    dynrules_output_rule
      (ptype, level,
       reaction_definition->reactant_definitions[i].rule_list, dynrules,
       dyncolors);
  }

  custom_printf(ptype, level, "\n");
}

void
reaction_definition_list_output(const output_type ptype, const int level,
                                reaction_definition_t *reaction_definition_list,
                                dynrules_t *dynrules, dyncolors_t *dyncolors)
{
  reaction_definition_t *current_reaction_definition;

  current_reaction_definition = reaction_definition_list;
  while (NULL != current_reaction_definition)
  {
    reaction_definition_output(ptype, level, current_reaction_definition,
                               dynrules, dyncolors);

    current_reaction_definition = current_reaction_definition->next;
  }
}

int reaction_definition_complete(reaction_definition_t *reaction_definition,
                                 dynrules_t *dynrules)
{
  int status;
  size_t used_rules_size;
  char *used_rules;

  used_rules_size = dynrules->count * sizeof(char);
  used_rules = malloc(used_rules_size);

  status = reaction_definition_complete_raw(reaction_definition, dynrules,
                                            used_rules_size, used_rules);

  free(used_rules);

  return status;
}

static int
reaction_definition_complete_raw(reaction_definition_t *reaction_definition,
                                 dynrules_t *dynrules,
                                 const size_t used_rules_size,
                                 char *used_rules)
{
  reactant_definition_t *current_reactant_definition;
  int i;

  current_reactant_definition = reaction_definition->reactant_definitions;
  for (i = 0; i < reaction_definition->reactant_definition_count;
       i++, current_reactant_definition++)
  {
    // Resolve rule name.
    current_reactant_definition->rule_list =
      dynrules->rules.a + current_reactant_definition->rule_index;
    if (NULL == current_reactant_definition->rule_list)
    {
      error_printf("Unable to find rule \"%d\" named in "
                   "reactant_definition \"%s\".\n",
                   current_reactant_definition->rule_index,
                   reaction_definition->name);
      return 0;
    }

    // Update *_counts
    memset(used_rules, 0, used_rules_size);

    reaction_definition->split_count =
      rule_get_max_primary_split_elements
        (current_reactant_definition->rule_list, used_rules);

    memset(used_rules, 0, used_rules_size);

    reaction_definition->bind_count =
      rule_get_max_bind_elements
        (current_reactant_definition->rule_list, used_rules);

    memset(used_rules, 0, used_rules_size);

    reaction_definition->self_bind_count =
      rule_get_max_primary_self_bind_elements
        (current_reactant_definition->rule_list, used_rules);

    memset(used_rules, 0, used_rules_size);

    reaction_definition->state_change_count =
      rule_get_max_state_change_elements
        (current_reactant_definition->rule_list, used_rules);

    memset(used_rules, 0, used_rules_size);

    reaction_definition->transient_bind_count =
      rule_get_max_transient_bind_elements
        (current_reactant_definition->rule_list, used_rules);

    memset(used_rules, 0, used_rules_size);

    reaction_definition->transient_self_bind_count =
      rule_get_max_primary_transient_self_bind_elements
        (current_reactant_definition->rule_list, used_rules);
  }

  // Create dynrule_aggregates to use during rule checking.
  reaction_definition->dynrule_aggregates =
    dynrule_aggregates_create
    (INITIAL_RULE_AGGREGATES,
     reaction_definition->split_count,
     reaction_definition->bind_count,
     reaction_definition->self_bind_count,
     reaction_definition->state_change_count,
     reaction_definition->transient_bind_count,
     reaction_definition->transient_self_bind_count,
     reaction_definition->product_count);

  return 1;
}

void
reaction_definition_find_max_rc(reaction_definition_t *reaction_definition_list,
                                double *pmax_uni, double *pmax_bi)
{
  reaction_definition_t *current_reaction_definition;
  double max_uni;
  double max_bi;

  max_uni = 0;
  max_bi = 0;

  current_reaction_definition = reaction_definition_list;
  while (NULL != current_reaction_definition)
  {
    if ((1 == current_reaction_definition->reactant_definition_count) &&
        (0 == current_reaction_definition->self_bind_count) &&
        (0 == current_reaction_definition->transient_self_bind_count))
    {
      // Unimolecular and uniparticle reaction.
      if (max_uni < current_reaction_definition->rate_const)
      {
        // Found larger unimolecular rate constant.
        max_uni = current_reaction_definition->rate_const;
      }
    }
    else
    {
      int self_bind_test;
      self_bind_test = current_reaction_definition->self_bind_count +
                       current_reaction_definition->transient_self_bind_count;

      if (((current_reaction_definition->reactant_definition_count > 2) ||
           (self_bind_test > 1)) ||
          ((current_reaction_definition->reactant_definition_count > 1) &&
           (self_bind_test != 0)))
      {
        error_printf("Do not support more than two reacting particles (%s).\n", current_reaction_definition->name);
        assert(0);
        exit(1);
      }

      // Bimolecular or biparticle reaction.
      if (max_bi < current_reaction_definition->rate_const)
      {
        // Found larger bimolecular or biparticle rate constant.
        max_bi = current_reaction_definition->rate_const;
      }
    }

    current_reaction_definition = current_reaction_definition->next;
  }

  *pmax_uni = max_uni;
  *pmax_bi = max_bi;
}

