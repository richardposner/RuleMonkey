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
#include <stdlib.h>
#include "dynps_reaction_containers.h"
#include "rule_aggregate.h"
#include "myrand.h"
#include "output.h"


void suppress_warnings39() {SUPPRESS_WARNINGS;}

// Local function declarations.
static boolean dynps_reaction_containers_resize
                 (dynps_reaction_containers_t *dynps_reaction_containers,
                  const size_t count);

static boolean dynps_reaction_containers_alloc
  (size_t *pnew_index, dynps_reaction_containers_t *dynps_reaction_containers);

static boolean dynps_reaction_containers_reaction_container_exists
  (dynps_reaction_containers_t *dynps_reaction_containers,
   reaction_definition_t  *reaction_definition,
   const int reactant_index1, int *aggregate1,
   const int reactant_index2, int *aggregate2);

void
dynps_reaction_containers_reset
  (dynps_reaction_containers_t *dynps_reaction_containers)
{
  dynps_reaction_containers->count = 0;
  dynps_reaction_containers->prob_sum = 0;
}

void
dynps_reaction_containers_output
  (const output_type ptype, const int level,
   dynps_reaction_containers_t *dynps_reaction_containers)
{
  int i;

  if (dynps_reaction_containers->count > 0)
  {
    CUSTOM_PRINTF(ptype, level, "prob_sum: %e\n",
                  dynps_reaction_containers->prob_sum);
  }

  for (i = 0; i < dynps_reaction_containers->count; i++)
  {
    CUSTOM_PRINTF(ptype, level, "%e: ",
                  dynps_reaction_containers->cumulative_probabilities.a[i]);
    PS_REACTION_CONTAINER_OUTPUT
      (ptype, level, &dynps_reaction_containers->ps_reaction_containers.a[i]);
  }
}

boolean
dynps_reaction_containers_reaction_container_add
  (dynps_reaction_containers_t *dynps_reaction_containers,
   reaction_definition_t *reaction_definition,
   const int reactant_index1, int *aggregate1,
   const int reactant_index2, int *aggregate2,
   const int full_aggregate_size)
{
  size_t new_index;
  ps_reaction_container_t *new_ps_reaction_container;
  double probability;
  boolean status;

  // Take into account any multipliers for probability.
  probability = reaction_definition->prob;

  if (0 == probability)
  {
    // Don't bother adding reaction if it can never occur.
    return TRUE;
  }

  // When adding A + B, make sure B + A isn't already added.
  if ((NULL != aggregate2) &&
      (1 == rule_aggregate_compare(aggregate1, aggregate2)) &&
      (TRUE == dynps_reaction_containers_reaction_container_exists
              (dynps_reaction_containers, reaction_definition,
               reactant_index2, aggregate2, reactant_index1, aggregate1)))
  {
    // Don't add A + B since B + A already exists.
    return TRUE;
  }
  
  status = dynps_reaction_containers_alloc
            (&new_index, dynps_reaction_containers);
  if (0 == status)
  {
    return FALSE;
  }

  if (NULL != aggregate1)
  {
    probability *= aggregate1[AGG_MULTIPLIER_OFFSET];
  }

  if (NULL != aggregate2)
  {
    probability *= aggregate2[AGG_MULTIPLIER_OFFSET];
  }

  new_ps_reaction_container =
    dynps_reaction_containers->ps_reaction_containers.a + new_index;

  ps_reaction_container_setup(full_aggregate_size, reaction_definition,
                              reactant_index1, aggregate1,
                              reactant_index2, aggregate2,
                              new_ps_reaction_container);

  // Calculate current probability sum.
  dynps_reaction_containers->prob_sum += probability;

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(3,
               "Adding reaction container '%s', prob:%le, agg1:%d, agg2:%d, "
               "prob_sum=%le.\n", reaction_definition->name,
               reaction_definition->prob,
               aggregate1[AGG_MULTIPLIER_OFFSET],
               aggregate2 ?
                 aggregate2[AGG_MULTIPLIER_OFFSET] : -1,
               dynps_reaction_containers->prob_sum);

#endif

  // Maintain list of cumulative probabilities so
  // we can choose between reactions.
  dynps_reaction_containers->cumulative_probabilities.a[new_index] =
    dynps_reaction_containers->prob_sum;

  return TRUE;
}

static boolean
dynps_reaction_containers_reaction_container_exists
  (dynps_reaction_containers_t *dynps_reaction_containers,
   reaction_definition_t *reaction_definition,
   const int reactant_index1, int *aggregate1,
   const int reactant_index2, int *aggregate2)
{ 
  size_t count;
  size_t sizeof_aggregate;
  size_t i;
  ps_reaction_container_t *current_container;

  count = dynps_reaction_containers->count;

  if (0 == count)
  {
    // Nothing to look for.
    return FALSE;
  }

  sizeof_aggregate = sizeof(int) *
    rule_aggregate_num_ints(
      reaction_definition->split_count,
      reaction_definition->bind_count,
      reaction_definition->self_bind_count,
      reaction_definition->state_change_count,
      reaction_definition->transient_bind_count,
      reaction_definition->transient_self_bind_count);

  for (i = 0,
       current_container = dynps_reaction_containers->ps_reaction_containers.a;
       i < count;
       i++, current_container++)
  {
    if ((reaction_definition->index ==
         current_container->reaction_definition->index) &&
        (reactant_index1 == current_container->reactant_index1) &&
        (reactant_index2 == current_container->reactant_index2) &&
        (0 == memcmp(aggregate1, current_container->aggregate1,
                     sizeof_aggregate)) &&
        (0 == memcmp(aggregate2, current_container->aggregate2,
                     sizeof_aggregate)))
    {
      // Full match found.
      return TRUE;
    }
  }

  // Match not found.
  return FALSE;
}

boolean
dynps_reaction_containers_attempt_reaction
  (dynps_reaction_containers_t *dynps_reaction_containers,
   ps_particle_t *particle1, ps_particle_t *particle2,
   world_t *world, dynps_species_t *dynps_species)
{
  size_t container_index;
  size_t count;
  double *current_cum_prob;
  boolean status;
  double my_rand;

  count = dynps_reaction_containers->count;
  assert(count != 0);
  container_index = 0;

  if (count > 1)
  {
    DYNPS_REACTION_CONTAINERS_OUTPUT(OUTPUT_DEBUG, 1,
                                     dynps_reaction_containers);
  }

  if (dynps_reaction_containers->prob_sum > 1)
  {
    specie_t *specie1;
    specie_t *specie2;
    specie1 = world->dynspecies->species.a + particle1->specie_index;
    if (NULL == particle2)
    {
      specie2 = NULL;
    }
    else
    {
      // Output more exhaustive error report.
      // #######*********######
    }
    DEBUG_PRINTF(0, "Probability too large: %e, count: %zd\n"
                 "Pass max_prob=>%e to simulate_ps.\n",
                 dynps_reaction_containers->prob_sum, count,
                 (1/dynps_reaction_containers->prob_sum) / count);
    error_printf("Probability too large: %e, count: %zd\n"
                 "Pass max_prob=>%e to simulate_ps.\n",
                 dynps_reaction_containers->prob_sum, count,
                 (1/dynps_reaction_containers->prob_sum) / count);
    dynps_reaction_containers_output(OUTPUT_ERROR, 0,
                                     dynps_reaction_containers);
    return FALSE;
  }

  status = TRUE;

  my_rand = MYRAND();

  for (current_cum_prob = dynps_reaction_containers->cumulative_probabilities.a;
       container_index < count;
       container_index++, current_cum_prob++)
  {
#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "  %lu :  current:%e, rand:%e\n", container_index,
           *current_cum_prob, my_rand);
#endif

    if (*current_cum_prob >= my_rand)
    {
      // Current index has been selected.

      status =
        ps_reaction_container_perform_reaction
        (dynps_reaction_containers->ps_reaction_containers.a + container_index,
         particle1, particle2, world, dynps_species);
      break;
    }
  }

  return status;
}

dynps_reaction_containers_t *
dynps_reaction_containers_create(const size_t init_count)
{
  dynps_reaction_containers_t *new_dynps_reaction_containers;
  int status;

  // All fields need to be initialized to zero anyways, so use calloc.
  new_dynps_reaction_containers =
    calloc(1, sizeof(dynps_reaction_containers_t));
  if (NULL == new_dynps_reaction_containers)
  {
    error_printf("Unable to allocate new dynps_reaction_containers.\n");
    return NULL;
  }

  status = dynarray_create2(init_count,
                            &new_dynps_reaction_containers->alloc_count,
                            sizeof(ps_reaction_container_t),
                            &new_dynps_reaction_containers->
                              ps_reaction_containers.v,
                            sizeof(double),
                            &new_dynps_reaction_containers->
                              cumulative_probabilities.v);
  if (0 == status)
  {
    // Error message output in dynarray_create.
    free(new_dynps_reaction_containers);
    return NULL;
  }

  return new_dynps_reaction_containers;
}

static boolean
dynps_reaction_containers_alloc
  (size_t *pnew_index, dynps_reaction_containers_t *dynps_reaction_containers)
{
  int status;

  status = dynps_reaction_containers_resize
             (dynps_reaction_containers, dynps_reaction_containers->count + 1);
  if (0 == status)
  {
    assert(0);
    return FALSE;
  }

  *pnew_index = dynps_reaction_containers->count++;

  return TRUE;
}

static boolean
dynps_reaction_containers_resize
  (dynps_reaction_containers_t *dynps_reaction_containers,
   const size_t count)
{
  int status;

  status = dynarray_resize2(count, &dynps_reaction_containers->alloc_count,
                            sizeof(ps_reaction_container_t),
                           &dynps_reaction_containers->
                             ps_reaction_containers.v,
                            sizeof(double),
                           &dynps_reaction_containers->
                             cumulative_probabilities.v);
  if (0 == status)
  {
    // Error message output in dynarray_resize.
    return FALSE;
  }

  return TRUE;
}

void
dynps_reaction_containers_destroy
  (dynps_reaction_containers_t *dynps_reaction_containers)
{
  free(dynps_reaction_containers->ps_reaction_containers.a);
  free(dynps_reaction_containers->cumulative_probabilities.a);
  free(dynps_reaction_containers);
}
