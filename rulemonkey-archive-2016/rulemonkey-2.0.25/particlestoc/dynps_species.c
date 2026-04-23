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
#include "dynps_species.h"
#include "dynarray.h"
#include "output.h"

void suppress_warnings36() {SUPPRESS_WARNINGS;}

// Local function declarations.
static boolean dynps_species_resize(dynps_species_t *dynps_species,
                                    const size_t count);
static boolean dynps_species_specie_create(int specie_index, cgraph_t *cg,
                                           dynps_species_t *dynps_species);
static size_t dynps_species_count_particles(dynps_species_t *dynps_species,
                                            dynspecies_t *dynspecies,
                                            dyncolors_t *dyncolors);

/**
** @brief Add specie to list.  Don't bother checking to see if equivalent
** specie is already in list.
**
**  @param cg Colored graph representation of specie
**  @param count Number of specie to add
**  @param update Whether or not specie should be updated when reaction occurs
**  @param last If adding multiple species, FALSE for all except TRUE for last
**  @param match If FALSE, must call ps_specie_check_all_reactants() and
**          report_definition_check_all() later
**  @param pspecie_id Will contain new specie_id if successfully
**  @param pnew Pointer to flag.  TRUE if new specie, otherwise will be FALSE
**  @param dynspecies Dynamic list of species
**  @param world Current state information of the simulated world
**  @param dynps_species Dynamic list of ps_species
*/
boolean
dynps_species_specie_add(cgraph_t *cg, const char *name,
                         const double pop_delta, boolean update,
                         boolean last, boolean match, int *pspecie_id,
                         boolean *pnew,
                         dynspecies_t *dynspecies, world_t *world,
                         dynps_species_t *dynps_species)
{
  int status;

#ifdef RULE_CHECK_DEBUG
      DEBUG_PRINTF(3, "In dynps_species_specie_add.\n");
#endif

  if (dynps_species->available.a != NULL)
  {
    // Previously freed specie is available for reuse.
    *pspecie_id = ilist_pop(&dynps_species->available.a);

#ifdef RULE_CHECK_DEBUG
      DEBUG_PRINTF(4, "Reusing specie %d.\n", *pspecie_id + 1);
#endif

    // Setup specie for use.
    status = specie_setup(*pspecie_id, cg, name, pop_delta, update, FALSE,
                          dynspecies->species.a + *pspecie_id);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }

    // Reset matching ps_specie.
    status = ps_specie_setup(cg, dynps_species->ps_species.a + *pspecie_id);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }
  }
  else
  {
    // Need to create new specie.
    status = dynspecies_specie_add(cg, name, pop_delta, update, last, match,
                                   FALSE, pspecie_id, pnew, dynspecies, world);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }

    // New specie was created, need to make matching ps_specie.
    status = dynps_species_specie_create(*pspecie_id, cg, dynps_species);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }
  }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(3, "Before rule check on new species:\n");
    SPECIE_OUTPUT(OUTPUT_DEBUG, 3,
                  dynspecies->species.a + *pspecie_id, world->dyncolors);
    PS_SPECIE_OUTPUT(OUTPUT_DEBUG, 3,
                     dynps_species->ps_species.a + *pspecie_id,
                     dynspecies->species.a + *pspecie_id);
#endif

  if (TRUE == match)
  {
    ps_specie_t *ps_specie;
    specie_t *specie;

    ps_specie = dynps_species->ps_species.a + *pspecie_id;
    specie = dynspecies->species.a + *pspecie_id;

    // Record which reactants new specie can participate as.
    status = ps_specie_check_all_reactants(ps_specie, specie, world);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }

    // Record which reports new specie should be included in.
    status = report_definition_check_all
               (specie, world->dynrules, world->report_definition_list, world);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(3, "Finished adding new ps_specie:\n");
    SPECIE_OUTPUT(OUTPUT_DEBUG, 3,
                  dynspecies->species.a + *pspecie_id, world->dyncolors);
    PS_SPECIE_OUTPUT(OUTPUT_DEBUG, 3,
                     dynps_species->ps_species.a + *pspecie_id,
                     dynspecies->species.a + *pspecie_id);
#endif
  }

  return TRUE;
}

/**
 ** @brief Move specie to unused list so that it can be reused later.
 **
 ** @param specie_id Specie to remove
 */
boolean
dynps_species_specie_remove(specie_t *specie, dynps_species_t *dynps_species)
{
  int status;

  ps_specie_reset(dynps_species->ps_species.a + specie->id);

  if (NULL != specie->name)
  {
    free(specie->name);
    specie->name = NULL;
  }

#ifdef RULE_CHECK_DEBUG
      DEBUG_PRINTF(4, "Removing specie %d.\n", specie->id + 1);
#endif

  status = ilist_add(specie->id, &dynps_species->available.a);
  if (status != 1)
  {
    return FALSE;
  }

  return TRUE;
}

/*
** @brief After initial species have been created, add all matching ps_species.
*/
boolean
dynps_species_specie_create_all(dynspecies_t *dynspecies,
                                dynps_species_t *dynps_species)
{
  boolean status;
  int i;
  int count;
  specie_t *species;

  count = dynspecies->count;
  species = dynspecies->species.a;
  for (i = dynps_species->count; i < count; i++)
  {
    specie_t *current_specie;

    current_specie = species + i;

    status = dynps_species_specie_create(i, current_specie->cg, dynps_species);
    if (FALSE == status)
    {
      return FALSE;
    }
  }

  return TRUE;
}

static boolean
dynps_species_specie_create(int specie_index, cgraph_t *cg,
                            dynps_species_t *dynps_species)
{
  boolean status;
  int ps_specie_index;

  // Need to create new specie.
  status = dynps_species_ps_specie_alloc(dynps_species);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

  ps_specie_index = dynps_species->count - 1;
  assert(ps_specie_index == specie_index);
  status = ps_specie_setup(cg, dynps_species->ps_species.a + ps_specie_index);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

  return TRUE;
}

boolean
dynps_species_check_all_reactants(dynps_species_t *dynps_species,
                                  dynspecies_t *dynspecies,
                                  world_t *world)
{
  int i;
  specie_t *current_specie;
  ps_specie_t *current_ps_specie;
  int count;
  int status;

  count = dynspecies->count;
  for (i = 0, current_ps_specie = dynps_species->ps_species.a,
              current_specie = dynspecies->species.a;
       i < count;
       i++, current_ps_specie++, current_specie++)
  {
    status = ps_specie_check_all_reactants(current_ps_specie, current_specie,
                                           world);
    if (FALSE == status)
    {
      return FALSE;
    }

    // Record which reports new specie should be included in.
    status = report_definition_check_all
               (current_specie, world->dynrules,
                world->report_definition_list, world);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }
  }

  return TRUE;
}

static size_t
dynps_species_count_particles(dynps_species_t *dynps_species,
                              dynspecies_t *dynspecies, dyncolors_t *dyncolors)
{
  size_t i;
  specie_t *current_specie;
  ps_specie_t *current_ps_specie;
  size_t count;
  size_t total_ps_particles;

  total_ps_particles = 0;

  count = dynspecies->count;
  for (i = 0, current_ps_specie = dynps_species->ps_species.a,
              current_specie = dynspecies->species.a;
       i < count;
       i++, current_ps_specie++, current_specie++)
  {
    size_t j;
    size_t particle_count;
    size_t active_particle_count;
#ifndef DS_NO_DOUBLECHECK
    size_t subtotal_ps_particles;
#endif

    subtotal_ps_particles = 0;
    particle_count = current_specie->cg->n;
    active_particle_count = 0;

    for (j = 0; j < particle_count; j++)
    {
      if (COLOR_TYPE_COMPONENT ==
          dyncolors_get_type_from_value(current_specie->cg->colors.a[j],
                                        dyncolors))
      {
        // Node is a component, so it is probably active.
        active_particle_count++;
      }
    }

    if (0 != current_specie->pop)
    {
#ifndef DS_NO_DOUBLECHECK
      subtotal_ps_particles = active_particle_count * current_specie->pop;
      assert((subtotal_ps_particles / current_specie->pop) ==
             active_particle_count);
      total_ps_particles += subtotal_ps_particles;
      assert(total_ps_particles >= subtotal_ps_particles);
#else
      total_ps_particles += active_particle_count * current_specie->pop;
#endif /* DS_NO_DOUBLECHECK */
    }
  }

  return total_ps_particles;
}

ps_particle_t *
dynps_species_setup_particles(dynspecies_t *dynspecies,
                              size_t *pps_particle_count,
                              dynps_species_t *dynps_species,
                              dyncolors_t *dyncolors)
{
  size_t i;
  specie_t *current_specie;
  ps_specie_t *current_ps_specie;
  size_t count;
  ps_particle_t *ps_particles;
  ps_particle_t *current_ps_particle;
  size_t ps_particle_count;

  ps_particle_count =
    dynps_species_count_particles(dynps_species, dynspecies, dyncolors);

  if (0 == ps_particle_count)
  {
    error_printf("No particles to create.\n");
    assert(FALSE);
    return NULL;
  }

  ps_particles = calloc(ps_particle_count, sizeof(ps_particle_t));
  if (NULL == ps_particles)
  {
    assert(FALSE);
    return NULL;
  }

  *pps_particle_count = ps_particle_count;

  current_ps_particle = ps_particles;

  count = dynspecies->count;
  for (i = 0, current_ps_specie = dynps_species->ps_species.a,
              current_specie = dynspecies->species.a;
       i < count;
       i++, current_ps_specie++, current_specie++)
  {
    double k;
    double current_specie_pop;
    size_t particle_count;

    current_specie_pop = current_specie->pop;
    particle_count = current_specie->cg->n;

    for (k = 0; k < current_specie_pop; k++)
    {
      size_t j;
      ps_particle_t *first_particle_in_specie;
      ps_particle_t *previous_particle_in_specie;

      // So that we can make linked list into ring.
      first_particle_in_specie = NULL;

      // So that we can create linked list.
      previous_particle_in_specie = NULL;

      for (j = 0; j < particle_count; j++)
      {
        if (COLOR_TYPE_COMPONENT ==
            dyncolors_get_type_from_value(current_specie->cg->colors.a[j],
                                          dyncolors))
        {
          ps_particle_setup
            (current_specie->id, j,
             current_ps_specie->dynps_reactant_containers_array[j],
             previous_particle_in_specie,
             current_ps_particle);

          if (NULL == first_particle_in_specie)
          {
            // Remeber first particle so that linked list can be made into ring.
            first_particle_in_specie = current_ps_particle;
          }

          previous_particle_in_specie = current_ps_particle;

          // Move on to next ps_particle.
          current_ps_particle++;
        }
      }

      if (first_particle_in_specie != NULL)
      {
        // Complete ring.
        first_particle_in_specie->next = previous_particle_in_specie;
      }
    }
  }

  return ps_particles;
}

void
dynps_species_output(const output_type ptype, const int level,
                     dynps_species_t *dynps_species,
                     dynspecies_t *dynspecies, dyncolors_t *dyncolors)
{
  int i;
  specie_t *current_specie;
  ps_specie_t *current_ps_specie;
  int count;

  count = dynspecies->count;
  for (i = 0, current_ps_specie = dynps_species->ps_species.a,
              current_specie = dynspecies->species.a;
       i < count;
       i++, current_ps_specie++, current_specie++)
  {
    SPECIE_OUTPUT(ptype, level, current_specie, dyncolors);
    PS_SPECIE_OUTPUT(ptype, level, current_ps_specie, current_specie);
  }
}

boolean
dynps_species_ps_specie_alloc(dynps_species_t *dynps_species)
{
  int status;

  status =
    dynps_species_resize(dynps_species, dynps_species->count + 1);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

  dynps_species->count++;

  return TRUE;
}

dynps_species_t *
dynps_species_create(const size_t init_count)
{
  dynps_species_t *new_dynps_species;
  int status;


  // All fields need to be initialized to zero anyways, so use calloc.
  new_dynps_species = calloc(1, sizeof(dynps_species_t));
  if (NULL == new_dynps_species)
  {
    error_printf("Unable to allocate new dynps_species.\n");
    return NULL;
  }

  // No unused species at the beginning.
  new_dynps_species->available.a = NULL;

  status = dynarray_create(init_count,
                           &new_dynps_species->alloc_count, sizeof(specie_t),
                           &new_dynps_species->ps_species.v);
  if (0 == status)
  {
    // Error message output in dynarray_create.
    free(new_dynps_species);
    return NULL;
  }

  return new_dynps_species;
}

static boolean
dynps_species_resize(dynps_species_t *dynps_species,
                  const size_t count)
{
  int status;

  status = dynarray_resize(count, &dynps_species->alloc_count,
                           sizeof(specie_t),
                           &dynps_species->ps_species.v);
  if (0 == status)
  {
    // Error message output in dynarray_resize.
    return FALSE;
  }

  return TRUE;
}

void
dynps_species_destroy(dynps_species_t *dynps_species)
{
  free(dynps_species->ps_species.a);
  free(dynps_species);
}
