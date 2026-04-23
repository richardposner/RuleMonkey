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

#include "ps_specie.h"
#include "ps_particle.h"
#include "naututil.h"
#include "constants.h"
#include "output.h"
#include "dynps_reaction_containers.h"
#include "rule_aggregate.h"

void suppress_warnings34() {SUPPRESS_WARNINGS;}

static boolean
ps_specie_check(ps_specie_t *ps_specie, specie_t *specie,
                dynrule_aggregates_t *dynrule_aggregates,
                reaction_definition_t *reaction_definition,
                int reactant_index, dynrules_t *dynrules,
                dyncolors_t *dyncolors);
static boolean
ps_specie_apply_dynrule_aggregates(ps_specie_t *ps_specie, specie_t *specie,
                               reaction_definition_t *reaction_definition,
                               int reactant_index,
                               dynrule_aggregates_t *dynrule_aggregates);

boolean
ps_specie_setup(cgraph_t *cg, ps_specie_t *new_ps_specie)
{
  size_t old_n;
  int i;

  if (new_ps_specie->dynps_reactant_containers_array != NULL)
  {
    assert(0);
    // Reusing pre-allocated object.
    old_n = new_ps_specie->n;
    new_ps_specie->n = cg->n;
    if (old_n < cg->n)
    {
      // Need more reactant containers.
      new_ps_specie->dynps_reactant_containers_array =
        realloc(new_ps_specie->dynps_reactant_containers_array,
                cg->n * sizeof(dynps_reactant_containers_t *));
      if (NULL == new_ps_specie->dynps_reactant_containers_array)
      {
        perror("Reallocating dynps_reactant_containers_array");
        assert(0);
        return FALSE;
      }
    }
  }
  else
  {
    // Creating new object.
    old_n = 0;
    new_ps_specie->n = cg->n;
    new_ps_specie->dynps_reactant_containers_array =
      calloc(cg->n, sizeof(dynps_reactant_containers_t *));
    if (NULL == new_ps_specie->dynps_reactant_containers_array)
    {
      perror("Allocating dynps_reactant_containers_array");
      assert(0);
      return FALSE;
    }
  }

  for (i = old_n; i < cg->n; i++)
  {
    new_ps_specie->dynps_reactant_containers_array[i] =
      dynps_reactant_containers_create(INITIAL_DYNREACTION_CONTAINERS);
    if (NULL == new_ps_specie->dynps_reactant_containers_array[i])
    {
      perror("Allocating dynps_reactant_containers_array element");
      assert(0);
      return FALSE;
    }
  }

  return TRUE;
}

/**
** @brief Reset data structure so that it can be reused.
**
** @param ps_specie data structure to reset.
*/
void
ps_specie_reset(ps_specie_t *ps_specie)
{
  int i;
  int n;
  dynps_reactant_containers_t **pdpsrc;

  n = ps_specie->n;

#if 0
  for (i = 0, pdpsrc = ps_specie->dynps_reactant_containers_array;
       i < n;
       i++, pdpsrc++)
  {
    dynps_reactant_containers_t *dpsrc;

    dpsrc = *pdpsrc;

    memset(dpsrc->ps_reactant_containers, 0,
           sizeof(ps_reactant_container_t) * dpsrc->count);
    dpsrc->count = 0;
  }
#else
  //if (NULL == ps_specie->dynps_reactant_containers_array)
  //{
   // return;
  //}

  for (i = 0, pdpsrc = ps_specie->dynps_reactant_containers_array;
       i < n;
       i++, pdpsrc++)
  {
    dynps_reactant_containers_destroy(*pdpsrc);
  }

  free(ps_specie->dynps_reactant_containers_array);
  ps_specie->dynps_reactant_containers_array = NULL;
  ps_specie->n = 0;
#endif
}

/* Given a specie, check all reactants for any to apply. */
boolean
ps_specie_check_all_reactants(ps_specie_t *ps_specie, specie_t *specie,
                              world_t *world)
{
  dynrule_aggregates_t *dynrule_aggregates;
  reaction_definition_t *current_reaction_definition;
  int status;

  /* Loop through every reactant in every reaction to find all matches */
  for (current_reaction_definition = world->reaction_definition_list;
       NULL != current_reaction_definition;
       current_reaction_definition = current_reaction_definition->next)
  {
    int reactant_index;
    int reactant_count;

    reactant_count =
     current_reaction_definition->reactant_definition_count;

    // Setup dynrule_aggregates to collect information about current reaction.
    dynrule_aggregates =
      dynrule_aggregates_create
        (INITIAL_RULE_AGGREGATES,
         current_reaction_definition->split_count,
         current_reaction_definition->bind_count,
         current_reaction_definition->self_bind_count,
         current_reaction_definition->state_change_count,
         current_reaction_definition->transient_bind_count,
         current_reaction_definition->transient_self_bind_count,
         current_reaction_definition->product_count);
    if (NULL == dynrule_aggregates)
    {
      error_printf("Error creating dynrule_aggregates.\n");
      return FALSE;
    }

    for (reactant_index = 0;
         reactant_index < reactant_count;
         reactant_index++)
    {
      dynrule_aggregates_reset_all(dynrule_aggregates);

      // Perform check.
      status = ps_specie_check(ps_specie, specie, dynrule_aggregates,
                               current_reaction_definition,
                               reactant_index, world->dynrules,
                               world->dyncolors);
      if (FALSE == status)
      {
        error_printf("Error calling ps_specie_check.\n");
        dynrule_aggregates_destroy(dynrule_aggregates);
        return FALSE;
      }
    }

    dynrule_aggregates_destroy(dynrule_aggregates);
  }

  return TRUE;
}

/* If the given specie matches the reactant definition, apply aggregate. */
static boolean
ps_specie_check(ps_specie_t *ps_specie, specie_t *specie,
                dynrule_aggregates_t *dynrule_aggregates,
                reaction_definition_t *reaction_definition,
                int reactant_index, dynrules_t *dynrules,
                dyncolors_t *dyncolors)
{
  int status;
  int count;

  status = specie_check(specie, dynrule_aggregates, reaction_definition,
                        reactant_index, dynrules, dyncolors, &count);
  if (0 == status)
  {
    error_printf("Error calling specie_check.\n");
    return FALSE;
  }

  if (count > 0)
  {
    // Specie was matched, so create matching particlestoc information.
    status = ps_specie_apply_dynrule_aggregates(ps_specie, specie,
                                            reaction_definition,
                                            reactant_index, dynrule_aggregates);
    if (FALSE == status)
    {
      // Current reactant doesn't do anything, just used for multiplicity.
      return TRUE;
    }
  }

  return TRUE;
}

static boolean
ps_specie_apply_dynrule_aggregates(ps_specie_t *ps_specie, specie_t *specie,
                               reaction_definition_t *reaction_definition,
                               int reactant_index,
                               dynrule_aggregates_t *dynrule_aggregates)
{
  int count;
  int i;
  int *current_aggregate;
  int pri;
  int sec;
  int id;
  int status;
  boolean bimolecular;

  count = dynrule_aggregates->count;
  assert(count > 0);

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "Applying:");
    DYNRULE_AGGREGATES_OUTPUT(OUTPUT_DEBUG, 5, dynrule_aggregates);
#endif
  for (i = 0,
       current_aggregate = dynrule_aggregates->rule_aggregates.a;
       i < count;
       i++, current_aggregate += dynrule_aggregates->rule_aggregate_num_ints)
  {
    size_t full_aggregate_size;

    full_aggregate_size = 
      dynrule_aggregates_sizeof_full_rule_aggregate(dynrule_aggregates);

    // Only note the first bind/split point.
    // This way probabalities are handled properly if more
    // than one connection is modified.
#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "first active of ");
    RULE_AGGREGATE_OUTPUT(OUTPUT_DEBUG, 5, current_aggregate);
#endif
    status = rule_aggregate_first_active(current_aggregate, &pri, &sec, &id,
					 &bimolecular);
    if (FALSE == status)
    {
      error_printf("\nError calling rule_aggregate_first_active for %s with:",
                   reaction_definition->name);
      RULE_AGGREGATE_OUTPUT(OUTPUT_ERROR, 0, current_aggregate);
      error_printf("\n");
      return FALSE;
    }
#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, ": pri:%d sec:%d", pri, sec);
    DEBUG_PRINTF(5, "\n");
#endif

    // Add primary bind/split point.
    status = dynps_reactant_containers_ps_reactant_container_add
               (reaction_definition, reactant_index, id,
                dynrule_aggregates->rule_aggregate_size, full_aggregate_size,
                current_aggregate,
                ps_specie->dynps_reactant_containers_array[pri]);
    if (FALSE == status)
    {
      DEBUG_PRINTF(2, "First active primary point not found.\n");
      return FALSE;
    }

    // Add secondary bind or self-bind  point if it exists.
    if ((-1 != sec) && (TRUE == bimolecular))
    {
      status = dynps_reactant_containers_ps_reactant_container_add
                 (reaction_definition, reactant_index, id,
                  dynrule_aggregates->rule_aggregate_size, full_aggregate_size,
                  current_aggregate,
                  ps_specie->dynps_reactant_containers_array[sec]);
      if (FALSE == status)
      {
        DEBUG_PRINTF(2, "First active secondary point not found.\n");
        return FALSE;
      }
    }
  }

  return TRUE;
}

boolean
ps_specie_attempt_reaction
  (ps_particle_t *particle1, ps_particle_t *particle2,
   ps_particle_t *ps_particles, world_t *world, dynps_species_t *dynps_species,
   dynps_reaction_containers_t *dynps_reaction_containers)
{
  int i;
  int j;
  int container1_count;
  int container2_count;
  dynps_reactant_containers_t *dynps_reactant_containers1;
  dynps_reactant_containers_t *dynps_reactant_containers2;
  ps_reactant_container_t *ps_reactant_container1;
  ps_reactant_container_t *ps_reactant_container2;
  int status;

#ifdef AGGRESSIVE_CHECK
  assert(particle1->specie_index < world->dynspecies->count);
  assert(particle1->specie_particle_index <
         world->dynspecies->species[particle1->specie_index].cg->n);
#endif

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(9, "Attempting reaction with the following particles:\n");
  PS_PARTICLE_OUTPUT(OUTPUT_DEBUG, 9, particle1, dynps_species, ps_particles);
  PS_PARTICLE_OUTPUT(OUTPUT_DEBUG, 9, particle2, dynps_species, ps_particles);
#endif

  dynps_reactant_containers1 = particle1->dynps_reactant_containers;
  container1_count = dynps_reactant_containers1->count;
  ps_reactant_container1 =
    dynps_reactant_containers1->ps_reactant_containers.a;

  if (NULL == particle2)
  {
    // Unimolecular reaction.

    for (i = 0; i < container1_count; i++, ps_reactant_container1++)
    {
      reaction_definition_t *reaction_definition;

      reaction_definition = ps_reactant_container1->reaction_definition;

      // If must have one reactant and
      // first particle is primary particle and
      // no self bindings and
      // no transient self bindings.
      if ((1 == reaction_definition->reactant_definition_count) &&
          ((ps_reactant_container1->primary_particle == -1) ||
           (ps_reactant_container1->primary_particle ==
           particle1->specie_particle_index)) &&
          (0 == reaction_definition->self_bind_count) &&
          (0 == reaction_definition->transient_self_bind_count))
      {
        // Reaction can be performed, so add to list.
        status = dynps_reaction_containers_reaction_container_add
          (dynps_reaction_containers,
           ps_reactant_container1->reaction_definition,
           ps_reactant_container1->reactant_index,
           ps_reactant_container1->aggregate,
           0,
           NULL,
           ps_reactant_container1->full_aggregate_size);

        if (FALSE == status)
        {
          return FALSE;
        }
      }
    }
  }
  else
  {
    // Bimolecular or biparticle reaction.
    dynps_reactant_containers2 = particle2->dynps_reactant_containers;
    container2_count = dynps_reactant_containers2->count;

    for (i = 0; i < container1_count; i++, ps_reactant_container1++)
    {
      for (j = 0, ps_reactant_container2 =
                  dynps_reactant_containers2->ps_reactant_containers.a;
           j < container2_count; j++, ps_reactant_container2++)
      {
        reaction_definition_t *reaction_definition1;
        reaction_definition_t *reaction_definition2;
        int reactant_index1;
        int reactant_index2;
        int reactant_pri1;
        int reactant_pri2;
        int p1_index;
        int p2_index;

        reaction_definition1 = ps_reactant_container1->reaction_definition;
        reaction_definition2 = ps_reactant_container2->reaction_definition;
        reactant_index1 = ps_reactant_container1->reactant_index;
        reactant_index2 = ps_reactant_container2->reactant_index;
        reactant_pri1 = ps_reactant_container1->primary_particle;
        reactant_pri2 = ps_reactant_container2->primary_particle;
        p1_index = particle1->specie_particle_index;
        p2_index = particle2->specie_particle_index;
        // If reaction definitions match and
        // primary particles match
        // then it might be a bimolecular or biparticle reaction.
        //
        // For bimolecular reaction:
        // Must have 2 reactants and
        // part of different reactants and
        // particles are not on same reactant.
        //
        // or
        //
        // For biparticle reaction:
        // primary particle is one of the chosen particles and
        // Must have 1 reactant and
        // must have 1 self bind point or transient self bind point and
        // particles must not map to same position on aggregate and
        // particles must belong to same aggregate and
        // both aggregates are identical.
        //
        // Then we have two components to perform a reaction.
        if ((reaction_definition1->index == reaction_definition2->index) &&
            (reactant_pri1 == reactant_pri2) &&
             // Bimolecular:
            (((2 == reaction_definition1->reactant_definition_count) &&
              (reactant_index1 != reactant_index2) &&
              (0 == ps_particle_check_same_entity(particle1, particle2))
             ) ||
             // Biparticle:
             ((1 == reaction_definition1->reactant_definition_count) &&
              ((reactant_pri1 == p1_index) || (reactant_pri1 == p2_index)) &&
              ((1 == reaction_definition1->self_bind_count) ||
               (1 == reaction_definition1->transient_self_bind_count)
              ) &&
              (p1_index != p2_index) &&
              (1 == ps_particle_check_same_entity(particle1, particle2)) &&
              (rule_aggregate_compare(ps_reactant_container1->aggregate,
                                      ps_reactant_container2->aggregate))
             )
            )
           )
        {
          // Reaction can be performed, so add to list.
#ifndef DS_NO_DOUBLECHECK
          // Both aggregates should have the same number of elements.
          assert(ps_reactant_container1->full_aggregate_size ==
                 ps_reactant_container2->full_aggregate_size);
#endif
          status = dynps_reaction_containers_reaction_container_add
            (dynps_reaction_containers,
             ps_reactant_container1->reaction_definition,
             ps_reactant_container1->reactant_index,
             ps_reactant_container1->aggregate,
             ps_reactant_container2->reactant_index,
             ps_reactant_container2->aggregate,
             ps_reactant_container1->full_aggregate_size);

          if (FALSE == status)
          {
            return FALSE;
          }
        }
      }
    }
  }

  // Attempt to perform one of the possible reactions.
  if (dynps_reaction_containers->count > 0)
  {
#ifdef RULE_CHECK_DEBUG
    DYNPS_REACTION_CONTAINERS_OUTPUT(OUTPUT_DEBUG, 5,
                                     dynps_reaction_containers);
#endif

    status =
      dynps_reaction_containers_attempt_reaction(dynps_reaction_containers,
                                                 particle1, particle2,
                                                 world, dynps_species);
    if (FALSE == status)
    {
      return FALSE;
    }
  }

  return TRUE;
}


void
ps_specie_output(const output_type ptype, const int level,
                 ps_specie_t *ps_specie, specie_t *specie)
{
  int i;
  dynps_reactant_containers_t **current_dynps_reactant_containers;

  CUSTOM_PRINTF(ptype, level, "ParticleStoc Reactant Containers:\n");
  current_dynps_reactant_containers = ps_specie->dynps_reactant_containers_array;
  if (NULL == current_dynps_reactant_containers)
  {
    CUSTOM_PRINTF(ptype, level, "  UNUSED SPECIE\n");
    return;
  }

  for (i = 0;
       i < specie->cg->n;
       i++, current_dynps_reactant_containers++)
  {
    CUSTOM_PRINTF(ptype, level, "  %d : ", i);
    if (NULL != *current_dynps_reactant_containers)
    {
     DYNPS_REACTANT_CONTAINERS_OUTPUT(ptype, level,
                                       *current_dynps_reactant_containers);
    }
    else
    {
      CUSTOM_PRINTF(ptype, level, "NULL\n");
    }
  }

  CUSTOM_PRINTF(ptype, level, "\n");
}
