/*
** Copyright 2008 The Translational Genomics Research Institute
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
** For each reaction_definition, we need to keep track of all the 
** specie permutations possible.  This is necessary for a network free
** Gillespie implementation.
*/
#ifndef RS_REACTION_DEFINITION_H
#define RS_REACTION_DEFINITION_H

typedef struct rs_reaction_definition_struct rs_reaction_definition_t;

#include "reaction_definition.h"
#include "constants.h"
#include "fsa.h"
#include "rs_specie.h"
#include "world.h"

#define RS_REACTION_DEFINITION_OUTPUT(ptype, level, rs_reaction_definition, \
                                      dynrules, dyncolors) \
   CHECKDEBUG(rs_reaction_definition_output, ptype, level, \
              rs_reaction_definition, dynrules, dyncolors)

struct rs_reaction_definition_struct
{
  // Pointer to reaction definition.
  reaction_definition_t *reaction_definition;

  // Used when applying rule.
  int *full_aggregate;

  /*
   * Computing the permutations is described below.
   *
   * Simply keep a list of species that interact as the first reactant and
   * a list of species that interact as the second reactant, and do not
   * pick the same specie instance for both first and second reactant.
   */
  // Doubly linked list of specie containers.
  fsa_t *rs_sc_fsa;


  // If 1, reaction is of type A(a) + A(a)
  // If 0, reaction is of type A(a) + A(b) or other unrelated.
  int same_components;

  // Sum of all the populations from each specie container for each type.
  // The following 4 variables only set when same_components == 1.
  // Total number of sites that just match first reactant only.
  double a_only_total;
  // Total number of sites that just match second reactant only.
  double b_only_total;
  // Total number of sites that match first and second reactant.
  double ab_both_total;

  // a_only_total + a_both_total.
  double a_total;
  // b_only_total + b_both_total.
  double b_total;

  // If same_components is 0, equal to the sum across each specie instance:
  // (# a sites that match 1st only) * (# b sites that match 2nd only) +
  // (# a sites that match 1st only) * (# b sites that match both) +
  // (# a sites that match both) * (# b sites that match 2nd only) +
  // (# a sites that match both) * (# b sites that match both).
  //
  // If same_components is 1, equal to the sum across each specie instance:
  // (# a sites that match 1st only) * (# b sites that match 2nd only) +
  // (# a sites that match 1st only) * (# b sites that both) +
  // (# a sites that match both) * (# b sites that match 2nd only) +
  // ((# a sites that match both) * ((# b sites that match both) - 1)) / 2.
  double total_extra;

  // For unimolecular, equal to a_total * reaction_definition->rate_constant.
  //
  // For bimolecular and same_components is 0, equal to
  // (a_total * b_total) - total_extra.
  //
  // For bimolecular and same_components is 1, equal to
  // (a_only_total * b_total) +
  // (a_both_total * b_only_total) +
  // ((a_both_total * (b_both_total - 1)) / 2) -
  // total_extra.
  double total_prob;
};


int rs_reaction_definition_perform_reaction
      (rs_reaction_definition_t *rs_reaction_definitions,
       size_t reaction_definition_count,
       fsa_t *rs_specie_fsa, world_t *world, double *ptotal_prob);

void
rs_reaction_definition_short_output(const output_type ptype, const int level,
                              rs_reaction_definition_t *rs_reaction_definition,
                              dynrules_t *dynrules, dyncolors_t *dyncolors);
int
rs_reaction_definition_init
  (rs_reaction_definition_t *rs_reaction_definition,
   reaction_definition_t *reaction_definition, dyncolors_t *dyncolors,
   dynrules_t *dynrules);

void rs_reaction_definition_output
  (const output_type ptype, const int level,
   rs_reaction_definition_t *rs_reaction_definition,
   dynrules_t *dynrules, dyncolors_t *dyncolors);

#endif /* RS_REACTION_DEFINITION_H */
