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

#ifndef REACTION_DEFINITION_H
#define REACTION_DEFINITION_H

typedef struct reaction_definition_struct reaction_definition_t;

#include "reactant_definition.h"
#include "llist.h"
#include "dynrule_aggregates.h"
#include "dynrules.h"
#include "dyncolors.h"
#include "output.h"

#define REACTION_FLAG_DIFFUSION 0x01

#define REACTION_DEFINITION_OUTPUT(ptype, level, ...) \
   CHECKDEBUG(reaction_definition_output, ptype, level, __VA_ARGS__)

#define REACTION_DEFINITION_LIST_OUTPUT(ptype, level, ...) \
   CHECKDEBUG(reaction_definition_list_output, ptype, level, __VA_ARGS__)

struct reaction_definition_struct
{
  int index;

  // For reporting.
  char *name;

  // Probability of reaction occuring.
  double rate_const;

  // Probability of reaction occuring.
  double prob;

  // Flags to change behavior of reaction.
  int flags;

  // Number of reactant definitions participating in reaction.
  int reactant_definition_count;

  // Number of seperate products created.
  int product_count;

  // Dynamically allocated array of pointers to
  // definitions of reactants that participate in reaction
  reactant_definition_t *reactant_definitions;

  // Number of splits that are registered.
  int split_count;

  // Number of binds that are registered.
  int bind_count;

  // Number of self binds that are registered.
  int self_bind_count;

  // Number of state changes that are registered.
  int state_change_count;

  // Number of transient binds that are registered.
  int transient_bind_count;

  // Number of transient self binds that are registered.
  int transient_self_bind_count;

  // Preallocated rule aggregate to use during rule checking.
  dynrule_aggregates_t *dynrule_aggregates;

  // List of species that need to be created.
  dynints_t *new_species;

  // For world linked list.
  reaction_definition_t *next;

  /* Variables below are responsible for tracking runtime state */

  // Dynamically allocated array of linked lists.
  // Each element in the array contains a list of reactants that
  // participates in the reaction.  An element will be NULL if no
  // species have matched the reactant rules yet.
  llist_t **reactant_llists;

  // Dynamically allocated array of linked lists.
  // Used to loop through all permutations of reactants.
  llist_t **reactant_work_llists;
  int *reactant_work_flags;
};

int reaction_definition_add(const char *name, const double prob,
                            const int flags, const int reactant_count,
                            const int product_count,
                            dynints_t *new_species,
                            reaction_definition_t **pnext);
void reaction_definition_output(const output_type ptype, const int level,
                                reaction_definition_t *reaction_definition,
                                dynrules_t *dynrules, dyncolors_t *dyncolors);
void reaction_definition_list_output
  (const output_type ptype, const int level,
   reaction_definition_t *reaction_definition_list, dynrules_t *dynrules,
   dyncolors_t *dyncolors);
int reaction_definition_complete(reaction_definition_t *reaction_definition,
                                 dynrules_t *dynrules);
void reaction_definition_find_max_rc
       (reaction_definition_t *reaction_definition_list,
        double *pmax_uni, double *pmax_multi);

#endif
