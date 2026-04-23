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


#include "dynps_species.h"
#include "ps_particle.h"
#include "constants.h"
#include "rule_aggregate.h"
#include "cgraph.h"
#include "output.h"
#include <assert.h>

void suppress_warnings37() {SUPPRESS_WARNINGS;}

static ps_particle_t *
ps_particle_renumber(ps_particle_t *beginning_particle1, const int *renum_a);

static boolean
ps_particle_cgraph_append_raw(ps_particle_t **pfull_particle,
                              cgraph_t *full_cgraph,
                              int *full_aggregate, ps_particle_t *particle1,
                              cgraph_t *cgraph1, int *aggregate1,
                              cgraph_t *final_cgraph,
                              int *renum_full, int *renum_1);
static boolean
ps_particle_add_new_species_raw(ps_particle_t *beginning_particle,
                                cgraph_t *full_cgraph,
                                cgraph_t *connected_cgraph,
                                cgraph_t *remaining_cgraph,
                                int *connected_renum, int *remaining_renum,
                                world_t *world, dynps_species_t *dynps_species);

void
ps_particle_setup(const size_t specie_index, const int specie_particle_index,
                  dynps_reactant_containers_t *dynps_reactant_containers,
                  ps_particle_t *next, ps_particle_t *new_ps_particle)
{
  new_ps_particle->specie_index = specie_index;
  new_ps_particle->specie_particle_index = specie_particle_index;
  new_ps_particle->dynps_reactant_containers = dynps_reactant_containers;
  new_ps_particle->next = next;
}

void
ps_particle_output_all(const output_type ptype, const int level,
                       ps_particle_t *ps_particles, int ps_particle_count,
                       dynps_species_t *dynps_species)
{
  int i;
  ps_particle_t *current_ps_particle;

  for (i = 0, current_ps_particle = ps_particles;
       i < ps_particle_count;
       i++, current_ps_particle++)
  {
    CUSTOM_PRINTF(ptype, level, "  %d : ", i);
    ps_particle_output(ptype, level, current_ps_particle, dynps_species,
                       ps_particles);
  }
}

void
ps_particle_output(const output_type ptype, const int level,
                   ps_particle_t *ps_particle, dynps_species_t *dynps_species,
                   ps_particle_t *ps_particles)
{
  ps_specie_t *current_ps_specie;
  dynps_reactant_containers_t *current_dynps_reactant_containers;

  if (NULL == ps_particle)
  {
    CUSTOM_PRINTF(ptype, level, "NULL\n");
    return;
  }

  CUSTOM_PRINTF(ptype, level, "%ld[%d]", ps_particle->specie_index + 1,
                ps_particle->specie_particle_index);
  if (NULL != ps_particles)
  {
    CUSTOM_PRINTF (ptype, level, " -> %ld",
                   (long)(ps_particle->next - ps_particles));
  }
  CUSTOM_PRINTF(ptype, level, " == ");
  current_ps_specie =
    dynps_species->ps_species.a + ps_particle->specie_index;
  current_dynps_reactant_containers =
    current_ps_specie->dynps_reactant_containers_array
    [ps_particle->specie_particle_index];

  if(current_dynps_reactant_containers !=
     ps_particle->dynps_reactant_containers)
  {
    CUSTOM_PRINTF (ptype, level, "<*mismatch*>");
    DYNPS_REACTANT_CONTAINERS_OUTPUT(ptype, level,
                                     current_dynps_reactant_containers);
    CUSTOM_PRINTF (ptype, level, "                   <* versus *>");
  }
  DYNPS_REACTANT_CONTAINERS_OUTPUT(ptype, level,
                                   ps_particle->dynps_reactant_containers);
}

void
ps_particle_output_specie(const output_type ptype, const int level,
                          ps_particle_t *first_ps_particle)
{
  ps_particle_t *ps_particle;

  ps_particle = first_ps_particle;

  CUSTOM_PRINTF(ptype, level, "%ld[%d]", ps_particle->specie_index + 1,
                ps_particle->specie_particle_index);

  ps_particle = ps_particle->next;
  while (ps_particle != first_ps_particle)
  {
    CUSTOM_PRINTF(ptype, level, ":");
    CUSTOM_PRINTF(ptype, level, "%ld[%d]", ps_particle->specie_index + 1,
           ps_particle->specie_particle_index);

    ps_particle = ps_particle->next;
  }

  CUSTOM_PRINTF(ptype, level, "\n");
}

boolean
ps_particle_check_same_entity(ps_particle_t *particle1,
                              ps_particle_t *particle2)
{
  ps_particle_t *current_particle;

  if (particle1 == particle2)
  {
    /* Same particle, so obviously same specie. */
    return TRUE;
  }

  current_particle = particle1->next;

  while(current_particle != particle1)
  {
    if (current_particle == particle2)
    {
      return TRUE;
    }

    current_particle = current_particle->next;
  }

  return FALSE;
}

boolean
ps_particle_cgraph_append(ps_particle_t **pfull_particle,
                          cgraph_t *full_cgraph,
                          int *full_aggregate,
                          ps_particle_t *particle1,
                          cgraph_t *cgraph1,
                          int *aggregate1,
                          cgraph_t *final_cgraph)
{
  int *renum_full; // Maps element numbers from a to c.
  int *renum_1; // Maps element numbers from b to c.
  int status;

  if (NULL == full_cgraph)
  {
    // Adding first reactant.
    renum_full = NULL;
    renum_1 = NULL;
  }
  else
  {
    renum_full = calloc((size_t)full_cgraph->n, sizeof(int));
    if (NULL == renum_full)
    {
      assert(0);
      return FALSE;
    }

    renum_1 = calloc((size_t)cgraph1->n, sizeof(int));
    if (NULL == renum_1)
    {
      assert(0);
      return FALSE;
    }
  }

  status = ps_particle_cgraph_append_raw(pfull_particle,
                                         full_cgraph,
                                         full_aggregate,
                                         particle1, cgraph1, aggregate1,
                                         final_cgraph,
                                         renum_full, renum_1);

  if (NULL != full_cgraph)
  {
    free(renum_full);
    free(renum_1);
  }

  return status;
}

static boolean
ps_particle_cgraph_append_raw(ps_particle_t **pfull_particle,
                              cgraph_t *full_cgraph,
                              int *full_aggregate, ps_particle_t *particle1,
                              cgraph_t *cgraph1, int *aggregate1,
                              cgraph_t *final_cgraph,
                              int *renum_full, int *renum_1)
{
  boolean status;

  status = rule_aggregate_cgraph_append_raw(full_cgraph, full_aggregate,
                                            cgraph1, aggregate1,
                                            final_cgraph,
                                            renum_full, renum_1);
  if (0 == status)
  {
    return FALSE;
  }

  if (NULL == *pfull_particle)
  {
    // Nothing in full_particle, so initialize with new particle.
    *pfull_particle = particle1;
  }
  else
  {
    ps_particle_t *last_full_particle;
    ps_particle_t *last_particle1;

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "About to combine the two particles:\n");
    PS_PARTICLE_OUTPUT_SPECIE(OUTPUT_DEBUG, 5, *pfull_particle);
    PS_PARTICLE_OUTPUT_SPECIE(OUTPUT_DEBUG, 5, particle1);
#endif

    // Renumber all particles.
    last_full_particle = ps_particle_renumber(*pfull_particle, renum_full);
    last_particle1 = ps_particle_renumber(particle1, renum_1);

    // Insert new particle chain into full_particle chain.
    last_full_particle->next = particle1;
    last_particle1->next = *pfull_particle;

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "Merged particles:\n");
    PS_PARTICLE_OUTPUT_SPECIE(OUTPUT_DEBUG, 5, *pfull_particle);
#endif
  }

  return TRUE;
}

static ps_particle_t *
ps_particle_renumber(ps_particle_t *beginning_particle1, const int *renum_a)
{
  ps_particle_t *previous_particle1;
  ps_particle_t *current_particle1;

  previous_particle1 = beginning_particle1;
  current_particle1 = beginning_particle1;

  for (;;)
  {
    current_particle1->specie_particle_index =
      renum_a[current_particle1->specie_particle_index];

    // Pointer to dynps_reactant_containers now invalid, will update once
    // Specie is known.
    current_particle1->dynps_reactant_containers = NULL;

    previous_particle1 = current_particle1;
    current_particle1 = current_particle1->next;
    if (current_particle1 == beginning_particle1)
    {
      // Renumbered all of the particles.
      // Return end of chain.
      return previous_particle1;
    }
  }
}

/*
** Note: After calling this function, the state of cgraph1 is indeterminate.
*/
boolean
ps_particle_add_new_species(ps_particle_t *beginning_particle,
                            cgraph_t *full_cgraph, world_t *world,
                            dynps_species_t *dynps_species)
{
  boolean status;
  int *connected_renum;
  int *remaining_renum;
  cgraph_t *connected_cgraph;
  cgraph_t *remaining_cgraph;

  connected_renum = calloc((size_t)full_cgraph->n, sizeof(int));
  remaining_renum = calloc((size_t)full_cgraph->n, sizeof(int));

  connected_cgraph = cgraph_create(INITIAL_CGRAPH_SIZE);
  remaining_cgraph = cgraph_create(INITIAL_CGRAPH_SIZE);

  status = ps_particle_add_new_species_raw(beginning_particle, full_cgraph,
                                           connected_cgraph, remaining_cgraph,
                                           connected_renum, remaining_renum,
                                           world, dynps_species);

  free(connected_renum);
  free(remaining_renum);

  cgraph_destroy(connected_cgraph);
  cgraph_destroy(remaining_cgraph);

  return status;
}

/**
 * @brief Reuse particles with the following new specie(s).
 */
static boolean
ps_particle_add_new_species_raw(ps_particle_t *beginning_particle,
                                cgraph_t *full_cgraph,
                                cgraph_t *connected_cgraph,
                                cgraph_t *remaining_cgraph,
                                int *connected_renum, int *remaining_renum,
                                world_t *world, dynps_species_t *dynps_species)
{
  boolean status;
  boolean new;
  ps_particle_t *current_particle;
  int *connected_lab;
  int *remaining_lab;
  int connected_specie_id;
  int remaining_specie_id;
  ps_particle_t *connected_particle_chain_begin;
  ps_particle_t *connected_particle_chain_end;
  ps_particle_t *remaining_particle_chain_begin;
  ps_particle_t *remaining_particle_chain_end;

  new = FALSE;
  connected_specie_id = 0;
  remaining_specie_id = 0;

  // Seperate fully connected graph from any other graphs that may exist.
  status = cgraph_seperate_raw(full_cgraph,
                               connected_cgraph, remaining_cgraph,
                               connected_renum, remaining_renum);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

  if (0 == connected_cgraph->n)
  {
    // Graph is empty.
    assert(0 == full_cgraph->n);
    assert(0 == remaining_cgraph->n);
    return TRUE;
  }

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(5, "Before calling dynps_species_specie_add\n");
  CGRAPH_OUTPUT(OUTPUT_DEBUG, 5, connected_cgraph, world->dyncolors);
#endif

  // Add first new specie.
  status = dynps_species_specie_add(connected_cgraph, NULL, 1, TRUE, TRUE,
                                    TRUE, &connected_specie_id, &new,
                                    world->dynspecies, world, dynps_species);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

#ifdef RULE_CHECK_DEBUG
  if (TRUE == new)
  {
    DEBUG_PRINTF(5, "Added new specie.\n");
  }
#endif

  if (remaining_cgraph->n != 0)
  {
    // Split was done, so need to add second specie.

    status = dynps_species_specie_add(remaining_cgraph, NULL, 1, TRUE, TRUE,
                                      TRUE, &remaining_specie_id, &new,
                                      world->dynspecies, world, dynps_species);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "Added second new specie.\n");
#endif

  }

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(5, "Beginning particle:\n");
  PS_PARTICLE_OUTPUT_SPECIE(OUTPUT_DEBUG, 5, beginning_particle);
  DEBUG_PRINTF(5, "connected_cgraph->n:%d\n", connected_cgraph->n);
  DEBUG_PRINTF(5, "remaining_cgraph->n:%d\n", remaining_cgraph->n);
#endif
  current_particle = beginning_particle;
  connected_lab = connected_cgraph->lab.a;
  remaining_lab = remaining_cgraph->lab.a;
  connected_particle_chain_begin = NULL;
  connected_particle_chain_end = NULL;
  remaining_particle_chain_begin = NULL;
  remaining_particle_chain_end = NULL;
  // Renumber all particles and update specie indexes.
  for (;;)
  {
    int old_specie_particle_index;
    int connected_specie_particle_index;
    int remaining_specie_particle_index;

    old_specie_particle_index = current_particle->specie_particle_index;
    connected_specie_particle_index =
      connected_renum[old_specie_particle_index];
    remaining_specie_particle_index =
      remaining_renum[old_specie_particle_index];

#ifdef AGGRESSIVE_CHECK
    assert(old_specie_particle_index < full_cgraph->n);
#endif

    if (connected_specie_particle_index != -1)
    {
#ifdef AGGRESSIVE_CHECK
      assert(-1 == remaining_specie_particle_index);
      assert(connected_specie_particle_index <
             connected_cgraph->n);
      assert(connected_lab[connected_specie_particle_index] <
             connected_cgraph->n);
#endif
      current_particle->specie_index = connected_specie_id;
      current_particle->specie_particle_index =
        connected_lab[connected_specie_particle_index];

      // Update pointer so to speed up reaction checking.
      current_particle->dynps_reactant_containers =
        dynps_species->ps_species.a[connected_specie_id].
        dynps_reactant_containers_array
          [connected_lab[connected_specie_particle_index]];
      assert(NULL != current_particle->dynps_reactant_containers);

      if (NULL == connected_particle_chain_begin)
      {
        connected_particle_chain_begin = current_particle;
      }
      else
      {
        connected_particle_chain_end->next = current_particle;
      }
      connected_particle_chain_end = current_particle;
    }
    else if (remaining_specie_particle_index != -1)
    {
#ifdef AGGRESSIVE_CHECK
      assert(-1 == connected_specie_particle_index);
      assert(remaining_specie_particle_index < remaining_cgraph->n);
      assert(remaining_lab[remaining_specie_particle_index] <
             remaining_cgraph->n);
#endif
      current_particle->specie_index = remaining_specie_id;
      current_particle->specie_particle_index =
        remaining_lab[remaining_specie_particle_index];

      // Update pointer so to speed up reaction checking.
      current_particle->dynps_reactant_containers =
        dynps_species->ps_species.a[remaining_specie_id].
        dynps_reactant_containers_array
          [remaining_lab[remaining_specie_particle_index]];
      assert(NULL != current_particle->dynps_reactant_containers);
      if (NULL == remaining_particle_chain_begin)
      {
        remaining_particle_chain_begin = current_particle;
      }
      else
      {
        remaining_particle_chain_end->next = current_particle;
      }
      remaining_particle_chain_end = current_particle;
    }
    else
    {
      error_printf("Error seperating disjoint cgraph.\n");
      assert(0);
      return FALSE;
    }

    current_particle = current_particle->next;
    if (current_particle == beginning_particle)
    {
      // Finished updating all modified particles.
      break;
    }
  }

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(5, "Updated particles:\n");
#endif

  // Complete circular chains of particles.
  if (NULL != connected_particle_chain_end)
  {
    connected_particle_chain_end->next = connected_particle_chain_begin;
#ifdef RULE_CHECK_DEBUG
    PS_PARTICLE_OUTPUT_SPECIE(OUTPUT_DEBUG, 5, connected_particle_chain_begin);
#endif
  }
  else
  {
    error_printf("Connected_particle_chain_end is not initialized.\n");
    assert(0);
    return FALSE;
  }

  if (NULL != remaining_particle_chain_end)
  {
    remaining_particle_chain_end->next = remaining_particle_chain_begin;
#ifdef RULE_CHECK_DEBUG
    PS_PARTICLE_OUTPUT_SPECIE(OUTPUT_DEBUG, 5, remaining_particle_chain_begin);
#endif
  }

  return TRUE;
}

