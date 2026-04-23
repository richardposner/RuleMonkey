/*
** Copyright 2008 Northern Arizona Univeristy and
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


typedef struct rs_specie_update_d_struct rs_specie_update_d_t;

#include "rs_specie.h"
#include "rs_specie_container.h"
#include <assert.h>
#include "output.h"
#include "world.h"
#include "specie.h"

static cgraph_t *connected_cgraph;
static cgraph_t *remaining_cgraph;
static cgraph_t *current_cgraph;

void suppress_warnings21() {SUPPRESS_WARNINGS;}

struct rs_specie_update_d_struct
{
  int delete;
  double pop_delta;
  double prob_delta;
};

static int rs_specie_add_new_species_raw
  (fsa_t *rs_specie_fsa,
   cgraph_t *full_cgraph,
   cgraph_t *connected_cgraph,
   cgraph_t *remaining_cgraph,
   cgraph_t *current_cgraph,
   world_t *world, double *pdelta_prob,
   rs_reaction_definition_t *rs_reaction_definitions);

static int rs_specie_update_cb(void *d, void *elem);

int
rs_specie_init(void *d, void *elem)
{
  int *pindex;
  int index;
  rs_specie_t *rs_specie;

  pindex = d;
  index = *pindex;
  rs_specie = elem;

  if (index >= 0)
  {
    // New rs_specie created, initilize index.
    rs_specie->specie_id = *pindex;
  }

  rs_specie->rs_specie_container_pointers_fsa =
    fsa_create(sizeof(rs_specie_container_t *), 8, NULL,
               rs_specie_container_pointer_clean_cb);
  if (NULL == rs_specie->rs_specie_container_pointers_fsa)
  {
    return 0;
  }

  return 1;
}

int
rs_specie_clean(void *d, void *elem)
{
  rs_specie_t *rs_specie;

  rs_specie = elem;

  // Delete rs_specie_containers.
  // Containers will be cleaned up automatically.
  fsa_destroy(rs_specie->rs_specie_container_pointers_fsa);

  // Important to set to NULL so that we know rs_specie is now unused.
  rs_specie->rs_specie_container_pointers_fsa = NULL;

  // rs_specie->index will be preserved so it can be reused later.

  return 1;
}

double
rs_specie_update(rs_specie_t *rs_specie, specie_t *specie, double pop_delta,
                 world_t *world, double *ptotal_prob)
{
  rs_specie_update_d_t rs_specie_update_d;
  int status;

  assert((specie->pop + pop_delta) >= 0);
  specie_update(specie, pop_delta, world);

  if ((specie->preserve == FALSE) && (specie->pop < 1))
  {
    // No species left, so delete all containers.
    rs_specie_update_d.delete = 1;
  }
  else
  {
    // Specie still populated, so preserve containers.
    rs_specie_update_d.delete = 0;
  }

  rs_specie_update_d.pop_delta = pop_delta;
  rs_specie_update_d.prob_delta = 0;

  // Update each reaction definition using specie containers.
  status = fsa_traverse(rs_specie_update_cb, &rs_specie_update_d,
                        rs_specie->rs_specie_container_pointers_fsa);
  if (0 == status)
  {
    return 0;
  }

  *ptotal_prob += rs_specie_update_d.prob_delta;

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(8, "rs_specie_update: pop_delta=%f; prob_delta=%f;"
               "new total_prob:%f\n", pop_delta,
               rs_specie_update_d.prob_delta, *ptotal_prob);
#endif

  return 1;
}

int
rs_specie_add_new_species(fsa_t *rs_specie_fsa, cgraph_t *full_cgraph,
                          world_t *world,
                          double *pdelta_prob,
                          rs_reaction_definition_t *rs_reaction_definitions)
{
  int status;

  if (NULL == connected_cgraph)
  {
    connected_cgraph = cgraph_create(INITIAL_CGRAPH_SIZE);
    remaining_cgraph = cgraph_create(INITIAL_CGRAPH_SIZE);
    current_cgraph = cgraph_create(INITIAL_CGRAPH_SIZE);
  }
  else
  {
    cgraph_reset(connected_cgraph);
    cgraph_reset(remaining_cgraph);
    cgraph_reset(current_cgraph);
  }

  status = rs_specie_add_new_species_raw(rs_specie_fsa, full_cgraph,
                                         connected_cgraph, remaining_cgraph,
                                         current_cgraph, world,
                                         pdelta_prob,
                                         rs_reaction_definitions);

  return status;
}

static int
rs_specie_add_new_species_raw
  (fsa_t *rs_specie_fsa,
   cgraph_t *full_cgraph,
   cgraph_t *connected_cgraph,
   cgraph_t *remaining_cgraph,
   cgraph_t *current_cgraph,
   world_t *world, double *pdelta_prob,
   rs_reaction_definition_t *rs_reaction_definitions)
{
  int status;
  cgraph_t *exchange_cgraph;
  dynspecies_t *dynspecies;
  specie_t *specie;
  rs_specie_t *rs_specie;
  boolean new;
  int product_count;

  status = cgraph_seperate(full_cgraph, connected_cgraph, remaining_cgraph);
  if (0 == status)
  {
    return 0;
  }

  assert(connected_cgraph->n > 0);

  dynspecies = world->dynspecies;

  product_count = 0;
  for (;;)
  {
    // Add new specie.
    product_count++;

#ifdef RULE_CHECK_DEBUG
    if (product_count > 1)
    {
      DEBUG_PRINTF(2, " + ");
    }
    CGRAPH_OUTPUT_BNGL(OUTPUT_DEBUG, 2, connected_cgraph, world->dyncolors);
#endif

    // Assume specie does not exist yet.
    new = 1;

#if 1 //#ifdef SPECIE_CHECK_DUP
    {
      // Search for existing specie.
      int i;
      specie_t *current_specie;
      rs_specie_t *current_rs_specie;
      int species_count;

      species_count = dynspecies->count;
      for (i = 0, current_specie = dynspecies->species.a;
           i < species_count;
           i++, current_specie++)
      {

        if (TRUE == cgraph_compare(connected_cgraph, current_specie->cg))
        {
          // Match found.
          // Update population.
          current_rs_specie =
            fsa_get_existing(current_specie->id, rs_specie_fsa);
          rs_specie_update(current_rs_specie, current_specie, 1.0, world,
                           pdelta_prob);

          new = 0;
          break;
        }
      }
    }
#endif

    if (1 == new)
    {
      rs_specie = fsa_get_new(rs_specie_fsa);
      if (NULL == rs_specie)
      {
        assert(0);
        return 0;
      }

      if (rs_specie->specie_id < dynspecies->count)
      {
        // Previously removed specie availabe to be reused.
        specie = dynspecies->species.a + rs_specie->specie_id;
        if (specie->name != NULL)
        {
          free(specie->name);
          specie->name = NULL;
        }
        status = specie_setup(rs_specie->specie_id, connected_cgraph,
                              "", 1, TRUE, FALSE, specie);
        if (0 == status)
        {
          return 0;
        }
      }
      else
      {
        int specie_id;

        // Need to create new specie.
        dynspecies_specie_add(connected_cgraph, "", 1, TRUE,
                              TRUE, TRUE, FALSE, &specie_id,
                              &new, dynspecies, world);
        if (FALSE == status)
        {
          return 0;
        }

        assert(specie_id == rs_specie->specie_id);

        specie = dynspecies->species.a + dynspecies->count - 1;
      }

      // Record which reactants new specie can participate in.
      status =
        rs_specie_check_all_reactants(rs_specie, rs_reaction_definitions,
                                      specie, world, pdelta_prob);
      if (FALSE == status)
      {
        return 0;
      }

      // Keep a note of which reports new specie should be included in.
      status =
        report_definition_check_all
          (specie, world->dynrules, world->report_definition_list, world);
      if (FALSE == status)
      {
        return FALSE;
      }
    }

    // Stop if no graph left.
    if (remaining_cgraph->n <= 0)
    {
      break;
    }

    // Seperate next graph.
    exchange_cgraph = current_cgraph;
    current_cgraph = remaining_cgraph;
    remaining_cgraph = exchange_cgraph;

    status = cgraph_seperate(current_cgraph,
                             connected_cgraph, remaining_cgraph);
    if (0 == status)
    {
      return 0;
    }
  }

  return 1;
}

void
rs_specie_output(const output_type ptype, const int level,
                 rs_specie_t *rs_specie, world_t *world)
{
  rs_specie_output_d_t rs_specie_output_d;

  rs_specie_output_d.ptype = ptype;
  rs_specie_output_d.level = level;
  rs_specie_output_d.world = world;

  rs_specie_output_cb(&rs_specie_output_d, rs_specie);
}

void
rs_specie_fsa_output(const output_type ptype, const int level,
                     fsa_t *rs_specie_fsa, world_t *world)
{
  rs_specie_output_d_t rs_specie_output_d;

  rs_specie_output_d.ptype = ptype;
  rs_specie_output_d.level = level;
  rs_specie_output_d.world = world;

  fsa_traverse(rs_specie_output_cb, &rs_specie_output_d, rs_specie_fsa);
}

int
rs_specie_output_cb(void *d, void *elem)
{
  rs_specie_output_d_t *rs_specie_output_d;
  rs_specie_t *rs_specie;
  world_t *world;

  rs_specie_output_d = d;
  rs_specie = elem;

  if (NULL == rs_specie->rs_specie_container_pointers_fsa)
  {
    // Not being used right now.
    return 1;
  }

  world = rs_specie_output_d->world;

  specie_output(rs_specie_output_d->ptype, rs_specie_output_d->level,
                world->dynspecies->species.a + rs_specie->specie_id,
                world->dyncolors);

  custom_printf(rs_specie_output_d->ptype, rs_specie_output_d->level,
                "Specie Containers:\n");
  fsa_traverse(rs_specie_container_pointer_output_cb, rs_specie_output_d,
               rs_specie->rs_specie_container_pointers_fsa);

  return 1;
}

void
rs_specie_remove(rs_specie_t *rs_specie, fsa_t *rs_specie_fsa)
{
  // Remove rs_specie so it can be reused later
  // (specie will be reused properly as well)
  fsa_remove(rs_specie, rs_specie_fsa);
}

static int
rs_specie_update_cb(void *d, void *elem)
{
  rs_specie_update_d_t *rs_specie_update_d;
  rs_specie_container_t **prs_specie_container;
  rs_specie_container_t *rs_specie_container;
  double prob_delta;

  rs_specie_update_d = d;
  prs_specie_container = elem;
  rs_specie_container = *prs_specie_container;

  prob_delta = rs_specie_container_update_pop
                 (rs_specie_container, rs_specie_update_d->pop_delta);

  // Keep running total of probability delta.
  rs_specie_update_d->prob_delta += prob_delta;

  return 1;
}

/* Given a specie, check all reactants for any to apply. */
int
rs_specie_check_all_reactants
  (rs_specie_t *rs_specie,
   rs_reaction_definition_t *rs_reaction_definitions, specie_t *specie,
   world_t *world, double *ptotal_prob)
{
  reaction_definition_t *current_reaction_definition;
  rs_reaction_definition_t *current_rs_reaction_definition;
  int status;
  double total_prob_delta;
  int count_array[2];
  dynrule_aggregates_t *dynrule_aggregates_array[2];
  set *set_array[2];
  int enable_rejection;

  total_prob_delta = 0;

  enable_rejection = world->flags & DS_FLAG_ENABLE_REJECTION;

  /* Loop through every reaction to find all matches */
  for (current_reaction_definition = world->reaction_definition_list,
       current_rs_reaction_definition = rs_reaction_definitions;
       NULL != current_reaction_definition;
       current_reaction_definition = current_reaction_definition->next,
       current_rs_reaction_definition++)
  {
    int reactant_count;
    int expected_product_count;
    rs_specie_container_t *rs_specie_container;
    double prob_delta;
    int i;

    reactant_count = current_reaction_definition->reactant_definition_count;
    expected_product_count = current_reaction_definition->product_count;
    count_array[0] = 0;
    count_array[1] = 0;
    dynrule_aggregates_array[0] = NULL;
    dynrule_aggregates_array[1] = NULL;

    for (i = 0; i < reactant_count; i++)
    {
      if (enable_rejection)
      {
        int node_count;
        int j;
        int rule_color;
        int *pcolor;

        rule_color = current_reaction_definition->reactant_definitions[i].
                     rule_list->color_value;
        // Simply find number of nodes that locally match.
        node_count = specie->cg->n;
        for (j = 0, pcolor = specie->cg->colors.a;
             j < node_count;
             j++, pcolor++)
        {
          if (*pcolor == rule_color)
          {
            // Node was locally matched.
            count_array[i]++;
            ADDELEMENT(set_array[i], j);
          }
        }

        if (count_array <= 0)
        {
          // No longer need set if nothing found..
          free(set_array[i]);
          set_array[i] = NULL;
        }
      }
      else
      {
        // Rejection free algorithm used.
        int product_count;

        dynrule_aggregates_array[i] =
          dynrule_aggregates_create
            (INITIAL_RULE_AGGREGATES,
             current_reaction_definition->split_count,
             current_reaction_definition->bind_count,
             current_reaction_definition->self_bind_count,
             current_reaction_definition->state_change_count,
             current_reaction_definition->transient_bind_count,
             current_reaction_definition->transient_self_bind_count,
             current_reaction_definition->product_count);
        if (NULL == dynrule_aggregates_array[i])
        {
          error_printf("Error creating dynrule_aggregates.\n");
          return 0;
        }

        // Perform check.
        status = specie_check(specie, dynrule_aggregates_array[i],
                              current_reaction_definition,
                              i, world->dynrules,
                              world->dyncolors, &count_array[i]);
        if (FALSE == status)
        {
          error_printf("Error calling specie_check.\n");
          return 0;
        }

        // Success. At this point, number of rule aggregates should be
        // the same as the match count.
        assert(count_array[i] == dynrule_aggregates_array[i]->count);

        // Check number of products produced (if split)
        product_count = dynrule_aggregates_array[i]->product_count;
        if (1 == reactant_count)
        {
          // Should always match.
          if ((0 != expected_product_count) &&
              (product_count != expected_product_count))
          {
            // Check failed.
            count_array[i] = 0;
          }
        }
        else if (product_count > 1)
        {
          // Don't support AB + C -> A + BC yet.
          error_printf("Don't support multiple bind/splits yet");
          return 0;
        }

        if (count_array[i] <= 0)
        {
          // No longer need dynrule_aggregates if nothing found.
          dynrule_aggregates_destroy(dynrule_aggregates_array[i]);
          dynrule_aggregates_array[i] = NULL;
        }
      }
    }

    if ((0 == count_array[0]) && (0 == count_array[1]))
    {
      // No reactants match.
      continue;
    }

    // Map reaction to specie and specie to reaction.
    rs_specie_container =
      rs_specie_container_create
        (current_rs_reaction_definition,
         specie, world->dyncolors, rs_specie, dynrule_aggregates_array,
         set_array, count_array);
    if (NULL == rs_specie_container)
    {
      error_printf("Error creating specie container\n");
      return 0;
    }

    // Initialize probabilities and permuations.
    prob_delta = rs_specie_container_update_pop
                   (rs_specie_container, specie->pop);

    // Keep running total of probabilities.
    total_prob_delta += prob_delta;
  }

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(8, "rs_specie_check_all_reactants: old_total_prob=%f; "
                  "total_prob_delta=%f\n", *ptotal_prob, total_prob_delta);
#endif
  *ptotal_prob += total_prob_delta;

  return 1;
} 

