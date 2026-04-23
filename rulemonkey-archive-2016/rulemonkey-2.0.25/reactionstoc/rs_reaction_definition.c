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


#include <assert.h>
#include "constants.h"
#include "rs_reaction_definition.h"
#include "rs_specie_container.h"
#include "output.h"
#include "myrand.h"
#include "rule_aggregate.h"
#include "dynps_species.h"
#include "dynrules.h"

void suppress_warnings40() {SUPPRESS_WARNINGS;}

typedef struct rs_specie_info_struct rs_specie_info_t;
struct rs_specie_info_struct
{
  // Which reactant we are looking for (input).
  int reactant_index;

  // Which one to pick (input and state tracking).
  double remaining_index;

  // Specie that has already been picked and should be ignored.
  // Set to -1 if no specie has been picked yet. (input)
  int ignore_specie_index;

  // Index of selected specie (output).
  int specie_index;
  // Selected rule aggregate (output).
  int *rule_aggregate;
  // Used when selecting 2nd reactant.
  dynrule_aggregates_t **dynrule_aggregates_array;
  int *count_array;

  // Used to check when same specie chosen for both reactants.
  int pop;
};

static int rs_reaction_definition_perform_reaction_raw
  (rs_specie_info_t *specie_info_a, rs_specie_info_t *specie_info_b,
   fsa_t *rs_sc_fsa, int *full_aggregate_arg, fsa_t *rs_specie_fsa,
   world_t *world, double *ptotal_prob,
   rs_reaction_definition_t *rs_reaction_definitions);
static int rs_reaction_definition_get_specie(void *d, void *elem);

int
rs_reaction_definition_init
  (rs_reaction_definition_t *rs_reaction_definition,
   reaction_definition_t *reaction_definition, dyncolors_t *dyncolors,
   dynrules_t *dynrules)
{
  size_t full_aggregate_size;

  rs_reaction_definition->reaction_definition = reaction_definition;

  full_aggregate_size = dynrule_aggregates_sizeof_full_rule_aggregate
                          (reaction_definition->dynrule_aggregates);

  rs_reaction_definition->full_aggregate = calloc(1, full_aggregate_size);
  if (NULL == rs_reaction_definition->full_aggregate)
  {
    return 0;
  }

  rs_reaction_definition->rs_sc_fsa = RS_SPECIE_CONTAINER_CREATE();
  if (NULL == rs_reaction_definition->rs_sc_fsa)
  {
    return 0;
  }

  if (2 == reaction_definition->reactant_definition_count)
  {
    reactant_definition_t *reactant_definitions;
    int status;

    if (2 == reaction_definition->product_count)
    {
      error_printf("Non-binding bimolecular reaction %s not supported\n",
                   reaction_definition->name);
      return 0;
    }

    reactant_definitions = reaction_definition->reactant_definitions;

#if 1
    DEBUG_PRINTF(2, "About to call dynrules_same_components\n");
    DYNRULES_OUTPUT_RULE(OUTPUT_DEBUG, 2, reactant_definitions[0].rule_list,
                         dynrules, dyncolors);
    DYNRULES_OUTPUT_RULE(OUTPUT_DEBUG, 2, reactant_definitions[1].rule_list,
                         dynrules, dyncolors);
                         
#endif
    status = dynrules_same_components(reactant_definitions[0].rule_list,
                                      reactant_definitions[1].rule_list,
                                      dyncolors);
    if (-1 == status)
    {
      return 0;
    }

    rs_reaction_definition->same_components = status;
  }
  else
  {
    rs_reaction_definition->same_components = 0;
  }

  return 1;
}

int
rs_reaction_definition_perform_reaction
  (rs_reaction_definition_t *rs_reaction_definitions,
   size_t reaction_definition_count, fsa_t *rs_specie_fsa,
   world_t *world, double *ptotal_prob)
{
  double current_prob_sum;
  double new_prob_sum;
  double chosen_prob;
  size_t i;
  rs_reaction_definition_t *rs_reaction_definition;
  int status;
  rs_specie_info_t specie_info_a;
  rs_specie_info_t specie_info_b;
  
  // Pick reaction to perform.
  chosen_prob = *ptotal_prob * MYRAND();

  new_prob_sum = 0;

#ifdef RULE_CHECK_DEBUG
      DEBUG_PRINTF(2, "Choosing reaction; total, chosen: %e, %e\n",
                   *ptotal_prob, chosen_prob);
#endif
  // Loop through every rs_reaction_definition until we find the right one.
  for (i = 0, rs_reaction_definition = rs_reaction_definitions;
       i < reaction_definition_count;
       i++, rs_reaction_definition++)
  {
    double current_total_prob;
    dynints_t *new_species;
             
    current_prob_sum = new_prob_sum;

    current_total_prob = rs_reaction_definition->total_prob;

    new_prob_sum = current_prob_sum + current_total_prob;

    if ((chosen_prob < new_prob_sum) && (current_total_prob > 0))
    {
      // Chosen reaction is within current reaction definition.
#ifdef RULE_CHECK_DEBUG

      DEBUG_PRINTF(2, "Performing reaction %d(%s)\n",
                   rs_reaction_definition->reaction_definition->index,
                   rs_reaction_definition->reaction_definition->name);
#endif

      // Choose first reactant.
      do
      {
        double new_b_total;

        specie_info_a.rule_aggregate = NULL;
        specie_info_b.rule_aggregate = NULL;
        specie_info_a.reactant_index = 0;
        specie_info_a.ignore_specie_index = -1;

        specie_info_a.remaining_index =
          MYRAND_INDEX(rs_reaction_definition->a_total);
        fsa_traverse(rs_reaction_definition_get_specie, &specie_info_a,
                     rs_reaction_definition->rs_sc_fsa);
        if (NULL == specie_info_a.rule_aggregate)
        {
          error_printf("Error getting 1st reactant.\n");
          rs_reaction_definition_output(OUTPUT_DEBUG, 0, 
                                        rs_reaction_definition,
                                        world->dynrules, world->dyncolors);
          assert(0);
        }

        if (rs_reaction_definition->
            reaction_definition->reactant_definition_count < 2)
        {
          // Only choosing one reactant, so we are done.
          break;
        }

        new_b_total =
          rs_reaction_definition->b_total - specie_info_a.count_array[1];
        if (new_b_total < 1)
        {
          // Specie picked for first reactant is only specie that can be
          // used as second reactant, so pick new a, ignoring specie just
          // picked.
          specie_info_a.rule_aggregate = NULL;
          specie_info_a.ignore_specie_index = specie_info_a.specie_index;

          // Choose first reactant again.
          specie_info_a.remaining_index =
            MYRAND_INDEX(rs_reaction_definition->a_total -
                         specie_info_a.count_array[0]);
          fsa_traverse(rs_reaction_definition_get_specie, &specie_info_a,
                       rs_reaction_definition->rs_sc_fsa);
          if (NULL == specie_info_a.rule_aggregate)
          {
            error_printf("Error getting 1st reactant.\n");
            rs_reaction_definition_output(OUTPUT_DEBUG, 0, 
                                          rs_reaction_definition,
                                          world->dynrules, world->dyncolors);
            assert(0);
          }
        }

        // Bimolecular reaction, choose second reactant.
        specie_info_b.reactant_index = 1;
        // Do not ignore specie already picked,
        // if same specie picked for 1st and 2nd reactant, must
        // choose both 1st and 2nd reactant again.
        specie_info_b.ignore_specie_index = -1;

        specie_info_b.remaining_index =
          MYRAND_INDEX(rs_reaction_definition->b_total);

        fsa_traverse(rs_reaction_definition_get_specie, &specie_info_b,
                     rs_reaction_definition->rs_sc_fsa);
        if (NULL == specie_info_b.rule_aggregate)
        {
          error_printf("Error getting 2nd reactant.\n");
          specie_info_b.reactant_index  = 1;
          rs_reaction_definition_output(OUTPUT_DEBUG, 0, 
                                        rs_reaction_definition,
                                        world->dynrules, world->dyncolors);
          assert(0);
        }

        // Keep picking a different pair of reactants if the same specie
        // was picked and there is only one of that specie.
      } while ((specie_info_a.specie_index == specie_info_b.specie_index) &&
               (specie_info_a.pop < 2));

#ifdef RULE_CHECK_DEBUG
      if (CUSTOM_OUTPUT_CHECK(2))
      {
        int j;
        rs_reaction_definition_t *rs_rd;

        DEBUG_PRINTF(2, "\n");
        for (j = i + 1, rs_rd = rs_reaction_definition + 1;
             j < reaction_definition_count;
             j++, rs_rd++)
        {
          current_prob_sum = new_prob_sum;

          current_total_prob = rs_rd->total_prob;

          new_prob_sum = current_prob_sum + current_total_prob;
        }
        DEBUG_PRINTF(2, "*** total_prob, final added prob: %e, %e\n",
                     *ptotal_prob, new_prob_sum);
      }
#endif


      // Reactants have been chosen, now perform reaction.
      status = rs_reaction_definition_perform_reaction_raw
                 (&specie_info_a, &specie_info_b,
                  rs_specie_fsa,
                  rs_reaction_definition->full_aggregate,
                  rs_specie_fsa,
                  world, ptotal_prob,
                  rs_reaction_definitions);
      if (0 == status)
      {
        return FALSE;
      }

      // Add any new species.
      new_species = rs_reaction_definition->reaction_definition->new_species;
      if (NULL != new_species)
      {
        int i;
        int count;
        int *pcurrent_specie_index;

        count = new_species->count;
        for (i = 0, pcurrent_specie_index = new_species->ints.a;
             i < new_species->count;
             i++, pcurrent_specie_index++)
        {
          rs_specie_t *rs_specie;
          specie_t *specie;

          specie = world->dynspecies->species.a + *pcurrent_specie_index;
          rs_specie = fsa_get_existing(*pcurrent_specie_index, rs_specie_fsa);
          status = rs_specie_update(rs_specie, specie, 1, world,
                                    ptotal_prob);
          if (0 == status)
          {
            error_printf("Error updating rs_specie 2\n");
            return 0;
          }
        }
      }

      // Reaction has been performed.
      return 1;
    }
  }

  // We dropped out without choosing anything.
  if (0 == new_prob_sum)
  {
    // *ptotal_prob was non-zero, but nothing left to do.
    debug_printf(1, "Dropped out of specie list without choosing anything\n"
                 "new_prob_sum: %e, total: %e\n", new_prob_sum, *ptotal_prob);

    *ptotal_prob = 0;
    return 1;
  }
  error_printf("Dropped out of specie list without choosing anything\n"
               "new_prob_sum: %e, total: %e\n", new_prob_sum, *ptotal_prob);
  assert(0);

  return FALSE;
}

static int
rs_reaction_definition_perform_reaction_raw
  (rs_specie_info_t *specie_info_a, rs_specie_info_t *specie_info_b,
   fsa_t *rs_sc_fsa, int *full_aggregate_arg, fsa_t *rs_specie_fsa,
   world_t *world, double *ptotal_prob,
   rs_reaction_definition_t *rs_reaction_definitions)
{
  specie_t *specie_a;
  rs_specie_t *rs_specie_a;
  specie_t *specie_b;
  rs_specie_t *rs_specie_b;
  cgraph_t *final_cgraph;
  int *full_aggregate;
  int status;
  double a_pop;
  double b_pop;

  specie_a = world->dynspecies->species.a + specie_info_a->specie_index;
  rs_specie_a = fsa_get_existing(specie_info_a->specie_index, rs_specie_fsa);
  a_pop = specie_a->pop;
  b_pop = 0;
  assert(a_pop > 0);

  // Apply reaction.
  if (NULL == specie_info_b->rule_aggregate)
  {
    // Unimolecular reaction.
    specie_b = NULL;
    rs_specie_b = NULL;
    final_cgraph = cgraph_copy(specie_a->cg);
    full_aggregate = specie_info_a->rule_aggregate;

    // Update population, but don't remove species until after reaction has
    // occured so that we can use aggregate.
    DEBUG_PRINTF(5, "Performing unimolecular reaction with specie %d(%f)\n",
                 specie_a->id + 1, a_pop);

#ifdef RULE_CHECK_DEBUG
    // Output reactant.
    CGRAPH_OUTPUT_BNGL(OUTPUT_DEBUG, 2, specie_a->cg, world->dyncolors);
#endif

    status = rs_specie_update(rs_specie_a, specie_a, -1, world, ptotal_prob);
    if (0 == status)
    {
      error_printf("Error updating rs_specie 1\n");
      return 0;
    }
  }
  else
  {
    // Bimolecular reaction.
    specie_b = world->dynspecies->species.a + specie_info_b->specie_index;
    rs_specie_b = fsa_get_existing(specie_info_b->specie_index, rs_specie_fsa);
    b_pop = specie_b->pop;
    assert(b_pop > 0);

    DEBUG_PRINTF(5, "Performing bimolecular reaction with species "
                 "%d(%f) %d(%f)\n",
                 specie_a->id + 1, a_pop,
                 specie_b->id + 1, b_pop);

#ifdef RULE_CHECK_DEBUG
    // Output reactants.
    CGRAPH_OUTPUT_BNGL(OUTPUT_DEBUG, 2, specie_a->cg, world->dyncolors);
    DEBUG_PRINTF(2, " + ");
    CGRAPH_OUTPUT_BNGL(OUTPUT_DEBUG, 2, specie_b->cg, world->dyncolors);
#endif

    // Create cgraph to store combined cgraph.
    final_cgraph = cgraph_create(INITIAL_CGRAPH_SIZE);
    if (NULL == final_cgraph)
    {
      assert(0);
      return 0;
    }

    // Use full aggregate as already allocated in rs_reaction_definition.
    full_aggregate = full_aggregate_arg;
    if (NULL == full_aggregate)
    {
      assert(0);
      return FALSE;
    }

    // Add delimiters to initialize full aggregate.
    rule_aggregate_init_full
      (full_aggregate, specie_info_a->rule_aggregate);

    // Merge first aggregate into full_aggregate.
    status = rule_aggregate_cgraph_append(NULL, full_aggregate,
                                          NULL, specie_info_a->rule_aggregate,
                                          NULL);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "Merged 1 to full_aggregate:\n");
    RULE_AGGREGATE_OUTPUT(OUTPUT_DEBUG, 5, specie_info_a->rule_aggregate);
    DEBUG_PRINTF(5, "\n");
    RULE_AGGREGATE_OUTPUT(OUTPUT_DEBUG, 5, full_aggregate);
    DEBUG_PRINTF(5, "\n");
#endif

    // Merge second aggregate into full aggregate and merge both cgraphs.
    status = rule_aggregate_cgraph_append(specie_a->cg,
                                          full_aggregate,
                                          specie_b->cg,
                                          specie_info_b->rule_aggregate,
                                          final_cgraph);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "Merged 2 to full_aggregate:\n");
    RULE_AGGREGATE_OUTPUT(OUTPUT_DEBUG, 5, specie_info_b->rule_aggregate);
    DEBUG_PRINTF(5, "\n");
    RULE_AGGREGATE_OUTPUT(OUTPUT_DEBUG, 5, full_aggregate);
    DEBUG_PRINTF(5, "\n");
#endif

    // Update population, but don't remove species until after reaction has
    // occured so that we can use aggregate.
    status = rs_specie_update(rs_specie_a, specie_a, -1, world,
                              ptotal_prob);
    if (0 == status)
    {
      error_printf("Error updating rs_specie 2\n");
      return 0;
    }
    status = rs_specie_update(rs_specie_b, specie_b, -1, world,
                              ptotal_prob);
    if (0 == status)
    {
      error_printf("Error in rs_specie_update 3\n");
      return 0;
    }
  }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "Before applying aggregate:\n");
    CGRAPH_OUTPUT(OUTPUT_DEBUG, 5, final_cgraph, world->dyncolors);
#endif

  // Perform reaction.
  status = rule_aggregate_cgraph_apply(final_cgraph, full_aggregate);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

#ifdef RULE_CHECK_DEBUG
  // Output result.
  DEBUG_PRINTF(2, " -> ");
#endif

  // We can remove species if needed now that reaction has finished.
  if ((specie_b != NULL) && (0 == b_pop) && (specie_a->id != specie_b->id))
  {
    rs_specie_remove(rs_specie_b, rs_specie_fsa);
  }
  if (0 == a_pop)
  {
    rs_specie_remove(rs_specie_a, rs_specie_fsa);
  }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "After applying aggregate:\n");
    CGRAPH_OUTPUT(OUTPUT_DEBUG, 5, final_cgraph, world->dyncolors);
#endif

  // Add product(s) to population, splitting cgraph if needed.
  status = rs_specie_add_new_species(rs_specie_fsa, final_cgraph, world,
                                     ptotal_prob, rs_reaction_definitions);
  if (FALSE == status)
  {
    error_printf("Error adding new species\n");
    return FALSE;
  }

  cgraph_destroy(final_cgraph);


  return TRUE;
}

/*
** Callback used to choose the specie that will participate in reaction.
*/
static int
rs_reaction_definition_get_specie(void *d, void *elem)
{
  rs_specie_info_t *rs_specie_info;
  rs_specie_container_t *rs_specie_container;
  int aggregate_index;
  double remaining_index;
  double factored_pop;
  dynrule_aggregates_t *dynrule_aggregates;
  int reactant_index;
  int specie_index;
  double pop;

  rs_specie_info = d;
  rs_specie_container = elem;

  pop = rs_specie_container->pop;

  if (0 == pop)
  {
    // Specie container was already deleted, so nothing to do.
    return 1;
  }

  remaining_index = rs_specie_info->remaining_index;
  reactant_index = rs_specie_info->reactant_index;

  specie_index = rs_specie_container->specie_index;
  factored_pop = rs_specie_container->factored_pop_array[reactant_index];
  if (0 == factored_pop)
  {
    // Nothing exists for current reactant.
    return 1;
  }

  dynrule_aggregates =
    rs_specie_container->dynrule_aggregates_array[reactant_index];

  if (rs_specie_info->ignore_specie_index == specie_index)
  {
    // Don't count reactant(s) from already chosen specie.
    factored_pop -= dynrule_aggregates->count;

    if (factored_pop < 1)
    {
      // Single member of current specie has already been selected, so ignore.
      return 1;
    }
  }

  if (remaining_index > factored_pop)
  {
    // Not in this specie.
    rs_specie_info->remaining_index -= factored_pop;
    return 1;
  }

  // Index is in this specie, find modulus to use as index to aggregate.
  aggregate_index =
    remaining_index - (dynrule_aggregates->count *
                       floor(remaining_index / dynrule_aggregates->count));
  assert(aggregate_index < dynrule_aggregates->count);

  rs_specie_info->rule_aggregate =
    dynrule_aggregates->rule_aggregates.a +
     (aggregate_index * dynrule_aggregates->rule_aggregate_num_ints);

  rs_specie_info->specie_index = rs_specie_container->specie_index;

  rs_specie_info->dynrule_aggregates_array =
    rs_specie_container->dynrule_aggregates_array;

  rs_specie_info->count_array = rs_specie_container->count_array;
  rs_specie_info->pop = pop;

  // Found specie, don't try any more.
  return 0;
}

void
rs_reaction_definition_short_output(const output_type ptype, const int level,
                              rs_reaction_definition_t *rs_reaction_definition,
                              dynrules_t *dynrules, dyncolors_t *dyncolors)
{
  reaction_definition_t *reaction_definition;

  if (!(rs_reaction_definition->a_total ||
        rs_reaction_definition->b_total))
  {
    // No reactants, so don't print anything.
    return;
  }

  reaction_definition = rs_reaction_definition->reaction_definition;

  custom_printf(ptype, level, "%d(%s)[%d]:",
                reaction_definition->index,
                reaction_definition->name,
                reaction_definition->reactant_definition_count);

  if (rs_reaction_definition->a_total)
  {
    custom_printf(ptype, level, " a:%.0f",
                  rs_reaction_definition->a_only_total);
  }
  if (rs_reaction_definition->b_total)
  {
    custom_printf(ptype, level, " b:%.0f",
                  rs_reaction_definition->b_only_total);
  }
  if (rs_reaction_definition->ab_both_total)
  {
    custom_printf(ptype, level, " both:%.0f",
                  rs_reaction_definition->ab_both_total);
  }
  custom_printf(ptype, level, "\n");
}

void
rs_reaction_definition_output(const output_type ptype, const int level,
                              rs_reaction_definition_t *rs_reaction_definition,
                              dynrules_t *dynrules, dyncolors_t *dyncolors)
{
  rs_specie_output_d_t rs_specie_output_d;

  rs_specie_output_d.ptype = ptype;
  rs_specie_output_d.level = level;

  reaction_definition_output(ptype, level,
                             rs_reaction_definition->reaction_definition,
                             dynrules, dyncolors);

  custom_printf(ptype, level, "(%s) only,total A:%f,%f, B:%f,%f; both:%f\n",
                rs_reaction_definition->same_components ? "same" : "diff",
                rs_reaction_definition->a_only_total,
                rs_reaction_definition->a_total,
                rs_reaction_definition->b_only_total,
                rs_reaction_definition->b_total,
                rs_reaction_definition->ab_both_total);

  custom_printf(ptype, level, "extra:%f\n",
                rs_reaction_definition->total_extra);

  custom_printf(ptype, level, "Total prob: %f\n",
                rs_reaction_definition->total_prob);

  custom_printf(ptype, level, "rs_specie_container list:\n");
  fsa_traverse(rs_specie_container_output_cb, &rs_specie_output_d,
               rs_reaction_definition->rs_sc_fsa);
}

