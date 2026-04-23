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
** Allows rs_reaction_definition keep track of associated species.
** This is necessary for a network free Gillespie implementation.
** Also used by species to update population counts in reaction definitions.
*/
#ifndef RS_SPECIE_CONTAINER_H
#define RS_SPECIE_CONTAINER_H

typedef struct rs_specie_container_struct rs_specie_container_t;

#define RS_SPECIE_CONTAINER_OUTPUT(ptype, level, rs_specie_container) \
   CHECKDEBUG(rs_specie_container_output, ptype, level, rs_specie_container)

#define RS_SPECIE_CONTAINER_CREATE() \
   FSA_CREATE(sizeof(rs_specie_container_t), 8, NULL, NULL)

#include "specie.h"
#include "rs_reaction_definition.h"
#include "constants.h"
#include "dynrule_aggregates.h"
#include "rs_specie.h"
#include "output.h"
#include "dyncolors.h"

struct rs_specie_container_struct
{
  // Index to specie.
  size_t specie_index;

  // Current population of specie, used to weight which specie participates
  // in reaction, so store here for quicker access.
  // Also need to factor in dynrule_aggregates->count below.
  double pop;

  // factored_pop_array[i] == pop * count_array[i];
  double factored_pop_array[2];

  // Sum of all the populations from each specie container for each type.
  // (# particles that only participate in 2nd reactant.)
  double a_only;
  // (# particles that only participate in 1st reactant.)
  double b_only;
  // (# particles that participate in both 1st and 2nd reactant.)
  // Note that when same_components = 1, a_both == b_both.
  double ab_both;

  // If rs_reaction_definition->same_components is 0,
  // equal to:
  // ((# a sites that match 1st only) * (# b sites that match 2nd only) +
  // a_only * b_only +
  // a_only * ab_both +
  // ab_both * b_only +
  // ab_both * ab_both.
  //
  // If rs_reaction_definition->same_components is 1,
  // equal to:
  // a_only * b_only +
  // a_only * ab_both +
  // ab_both * b_only +
  // (ab_both * (ab_both - 1) / 2.
  double extra;

  // List of rule aggregates that fully defines each reaction that can occur.
  // Dynrule_aggregate(s) and array needs to be cleaned up when container
  // is deleted.
  // Is empty if partial check was done.
  dynrule_aggregates_t *dynrule_aggregates_array[2];

  // If bimolecular: array of elemenets for each aggregate
  // which contains an array for each RHS specie that is bound to
  // which contains an array for each aggregate in RHS specie
  // which contains an index to the specie(s) created.
  // Only the LHS reactant stores this information.
  //dynarray_t *agg_rspecie_ragg_products_array;

  // If unimolecular: array of elements for each aggregate
  // which contains index(es) to the specie(s) created.
  //dynints_t *agg_products_array;

  // List of sets that identify the nodes that were partially matched.
  // Is empty if full check was done.
  // Not used right now.
  //set **set_array;

  // Is set to the number of nodes matched. If dynrule_aggregates_array
  // is not NULL, count_array[i] is same as dynrule_aggregates_array[i]->count.
  // count_array[0] == a_total == a_only + ab_both
  // count_array[1] == b_total == b_only + ab_both
  int count_array[2];

  // Pointer to rs_reaction_definition, used when updating population and
  // permutation counts (do not delete when container is deleted).
  rs_reaction_definition_t *rs_reaction_definition;
};

rs_specie_container_t *rs_specie_container_create
  (rs_reaction_definition_t *rs_reaction_definition,
   specie_t *specie, dyncolors_t *dyncolors, rs_specie_t *rs_specie,
   dynrule_aggregates_t **dynrule_aggregates_array, set **set_array,
   int *count_array);
int rs_specie_container_pointer_clean_cb(void *d, void *elem);
int rs_specie_container_pointer_output_cb(void *d, void *elem);
int rs_specie_container_output_cb(void *d, void *elem);
void rs_specie_container_output(output_type  ptype, int level,
                                rs_specie_container_t *rs_specie_container);
double rs_specie_container_update_pop
  (rs_specie_container_t *rs_specie_container, double pop_delta);

#endif /* RS_SPECIE_CONTAINER_H */

