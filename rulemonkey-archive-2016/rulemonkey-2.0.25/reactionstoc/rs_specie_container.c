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


#include <assert.h>
#include "dynrule_aggregates.h"
#include "rs_specie_container.h"
#include "constants.h"

void suppress_warnings17() {SUPPRESS_WARNINGS;}

int
rs_specie_container_pointer_clean_cb(void *d, void *elem)
{
  rs_specie_container_t **prs_specie_container;
  rs_specie_container_t *rs_specie_container;
  dynrule_aggregates_t **dynrule_aggregates_array;

  prs_specie_container = elem;
  rs_specie_container = *prs_specie_container;
  dynrule_aggregates_array = rs_specie_container->dynrule_aggregates_array;

  fsa_remove(rs_specie_container,
             rs_specie_container->rs_reaction_definition->rs_sc_fsa);

  if (NULL != dynrule_aggregates_array[0])
  {
    dynrule_aggregates_destroy(dynrule_aggregates_array[0]);
    dynrule_aggregates_array[0] = NULL;
  }
  if (NULL != dynrule_aggregates_array[1])
  {
    dynrule_aggregates_destroy(dynrule_aggregates_array[1]);
    dynrule_aggregates_array[0] = NULL;
  }

  rs_specie_container->pop = 0;
  rs_specie_container->factored_pop_array[0] = 0;
  rs_specie_container->factored_pop_array[1] = 0;

  rs_specie_container->rs_reaction_definition = NULL;

  return 1;
}

int
rs_specie_container_pointer_output_cb(void *d, void *elem)
{   
  rs_specie_output_d_t *rs_specie_output_d;
  rs_specie_container_t **prs_specie_container;
  rs_specie_container_t *rs_specie_container;
  output_type ptype;
  int level;
               
  rs_specie_output_d = d;
  prs_specie_container = elem;
  rs_specie_container = *prs_specie_container;

  ptype = rs_specie_output_d->ptype;
  level = rs_specie_output_d->level;

  rs_specie_container_output(ptype, level, rs_specie_container);

  return 1;
} 
  
int
rs_specie_container_output_cb(void *d, void *elem)
{
  rs_specie_output_d_t *rs_specie_output_d;
  rs_specie_container_t *rs_specie_container;
  output_type ptype;
  int level;

  rs_specie_output_d = d;
  rs_specie_container = elem;

  ptype = rs_specie_output_d->ptype;
  level = rs_specie_output_d->level;

  rs_specie_container_output(ptype, level, rs_specie_container);

  return 1;
}

void
rs_specie_container_output(output_type ptype, int level,
                           rs_specie_container_t *rs_specie_container)
{
  if (NULL == rs_specie_container->rs_reaction_definition)
  {
    // Specie container not in use right now, so nothing to output.
    return;
  }

  custom_printf(ptype, level, "%s %ld(%f[%f;%f]*[%d;%d]){%f;%f;%f}%f\n",
                rs_specie_container->rs_reaction_definition->
                  reaction_definition->name,
                rs_specie_container->specie_index + 1,
                rs_specie_container->pop,
                rs_specie_container->factored_pop_array[0],
                rs_specie_container->factored_pop_array[1],
                rs_specie_container->count_array[0],
                rs_specie_container->count_array[1],
                rs_specie_container->a_only,
                rs_specie_container->b_only,
                rs_specie_container->ab_both,
                rs_specie_container->extra);

  return;
}

/**
** Creates a specie container and initializes everything inside it
**
** Note: does not update rs_reaction_definition permutations or probabilities.
**       Call rs_specie_container_update_pop to do that.
*/
rs_specie_container_t *
rs_specie_container_create
  (rs_reaction_definition_t *rs_reaction_definition,
   specie_t *specie, dyncolors_t *dyncolors, rs_specie_t *rs_specie,
   dynrule_aggregates_t **dynrule_aggregates_array, set **set_array,
   int *count_array)
{
  rs_specie_container_t *rs_specie_container;
  rs_specie_container_t **rs_sc_p;
  // Store new rs_specie_container into rs_reaction_definition.
  rs_specie_container = fsa_get_new(rs_reaction_definition->rs_sc_fsa);

  rs_specie_container->specie_index = specie->id;
  // pop will be updated in rs_specie_container_update_pop.
  rs_specie_container->rs_reaction_definition = rs_reaction_definition;
  memcpy(rs_specie_container->dynrule_aggregates_array,
         dynrule_aggregates_array, 2 * sizeof(dynrule_aggregates_t *));
  //rs_specie_container->set_array = set_array;
  memcpy(rs_specie_container->count_array, count_array, 2 * sizeof(int));

  if (2 == rs_reaction_definition->
      reaction_definition->reactant_definition_count)
  {
    double  a_only;
    double  b_only;
    double  ab_both;

    if ((NULL != dynrule_aggregates_array[0]) &&
        (NULL != dynrule_aggregates_array[1]))
    {
      // Specie populated 1st and 2nd reactants, find number of
      // shared components so that probability can be computed.
      dynrule_aggregates_count_shared_components
        (dynrule_aggregates_array[0], dynrule_aggregates_array[1],
         specie, dyncolors, &a_only, &b_only, &ab_both);

      assert(count_array[0] == a_only + ab_both);
      assert(count_array[1] == b_only + ab_both);
    }
    else
    {
      // Specie only populated 1st or 2nd reactant.
      a_only = count_array[0];
      b_only = count_array[1];
      ab_both = 0;
    }

    rs_specie_container->extra =
      (a_only * b_only) + (a_only * ab_both) + (ab_both * b_only);
    if (0 == rs_reaction_definition->same_components)
    {
      rs_specie_container->extra += ab_both * ab_both;
    }
    else
    {
      rs_specie_container->extra += floor((ab_both * (ab_both - 1)) / 2);
    }

    rs_specie_container->a_only = a_only;
    rs_specie_container->b_only = b_only;
    rs_specie_container->ab_both = ab_both;
  }
  else
  {
    rs_specie_container->a_only = dynrule_aggregates_array[0]->count;
    rs_specie_container->b_only = 0;
    rs_specie_container->ab_both = 0;

    rs_specie_container->extra = 0;
  }

  // Store rs_specie_container in specie so that populations can be updated.
  rs_sc_p = fsa_get_new(rs_specie->rs_specie_container_pointers_fsa);
  *rs_sc_p = rs_specie_container;

  return rs_specie_container;
}

/**
** Updates population totals and permutations
*/
double
rs_specie_container_update_pop
  (rs_specie_container_t *rs_specie_container, double pop_delta)
{
  rs_reaction_definition_t *rs_reaction_definition;
  dynrule_aggregates_t **dynrule_aggregates_array;
  double rate_constant;
  double new_prob;
  double prob_delta;
  double new_pop;
  int reactant_count;
  int *count_array;
  double a_total;
  double total_extra;

  count_array = rs_specie_container->count_array;
  rs_reaction_definition = rs_specie_container->rs_reaction_definition;
  dynrule_aggregates_array = rs_specie_container->dynrule_aggregates_array;
  rate_constant = rs_reaction_definition->reaction_definition->rate_const;
  reactant_count =
    rs_reaction_definition->reaction_definition->reactant_definition_count;

  new_pop = rs_specie_container->pop + pop_delta;
  new_prob = 0;

  a_total = rs_reaction_definition->a_total + (count_array[0] * pop_delta);
  rs_specie_container->factored_pop_array[0] = new_pop * count_array[0];

  if (1 == rs_reaction_definition->reaction_definition->
      reactant_definition_count)
  {
    // Unimolecular.
    new_prob = a_total * rate_constant;
  }
  else
  {
    double b_total;

    b_total = rs_reaction_definition->b_total + (count_array[1] * pop_delta);
    total_extra = rs_reaction_definition->total_extra +
                  (rs_specie_container->extra * pop_delta);
    rs_specie_container->factored_pop_array[1] = new_pop * count_array[1];

    if (0 == rs_reaction_definition->same_components)
    {
      // Bimolecular, binding different components.
      new_prob = ((a_total * b_total) - total_extra) * rate_constant;
    }
    else
    {
      double a_only_total;
      double b_only_total;
      double ab_both_total;

      // Bimolecular, binding same components.
      a_only_total = rs_reaction_definition->a_only_total +
                     (rs_specie_container->a_only * pop_delta);

      b_only_total = rs_reaction_definition->b_only_total +
                     (rs_specie_container->b_only * pop_delta);
      ab_both_total = rs_reaction_definition->ab_both_total +
                     (rs_specie_container->ab_both * pop_delta);

      new_prob = ((a_only_total * b_total) +
                  (ab_both_total * b_only_total) +
                  ((ab_both_total * (ab_both_total - 1)) / 2) -
                  total_extra) *
                 rate_constant;

      rs_reaction_definition->a_only_total = a_only_total;
      rs_reaction_definition->b_only_total = b_only_total;
      rs_reaction_definition->ab_both_total = ab_both_total;
    }

    rs_reaction_definition->b_total = b_total;
    rs_reaction_definition->total_extra = total_extra;
  }

  rs_reaction_definition->a_total = a_total;

  rs_specie_container->pop = new_pop;

  prob_delta = new_prob - rs_reaction_definition->total_prob;

  rs_reaction_definition->total_prob = new_prob;

  return prob_delta;
}

