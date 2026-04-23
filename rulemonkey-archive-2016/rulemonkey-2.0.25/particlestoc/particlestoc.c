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
#include "particlestoc.h"
#include "constants.h"
#include "ps_reaction_definition.h"
#include "ps_particle.h"
#include "myrand.h"
#include "dynps_species.h"
#include "world.h"
#include "output.h"
#include "math.h"

void suppress_warnings35() {SUPPRESS_WARNINGS;}

// Local function declarations.
static double particlestoc_get_pseudo_particle_count
  (reaction_definition_t *reaction_definition_list,
   const size_t particle_count, const double volume,
   const double pseudo_factor,
   const double pseudo_count, const double avogadro);

/* \brief perform all initial setup needed before simulation can be run
**
** \param world State information
** \param dynps_species State information
** \param pps_particles array of ps_particles
** \param pparticle_count Pointer to current particle count (may be modified)
*/
boolean
particlestoc_complete(world_t *world, dynps_species_t *dynps_species,
                      ps_particle_t **pps_particles, size_t *pparticle_count)
{
  boolean status;

  status = world_complete(world);
  if (FALSE == status)
  {
    error_printf("Error setting up world.\n");
    return FALSE;
  }

  status = dynps_species_specie_create_all(world->dynspecies, dynps_species);
  if (FALSE == status)
  {
    error_printf("Error creating ps species.\n");
    return FALSE;
  }

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(3, "About to check all reactants.\n");
  DEBUG_PRINTF(3, "All reaction definitions:\n");
  REACTION_DEFINITION_LIST_OUTPUT(OUTPUT_DEBUG, 3,
                                  world->reaction_definition_list,
                                  world->dynrules, world->dyncolors);
  DEBUG_PRINTF(3, "Done with reaction definitions.\n");
  DEBUG_PRINTF(3, "All report definitions:\n");
  REPORT_DEFINITION_LIST_OUTPUT(OUTPUT_DEBUG, 3, world);
  DEBUG_PRINTF(3, "Done with report definitions.\n");
#endif

  status = dynps_species_check_all_reactants(dynps_species, world->dynspecies,
                                             world);
  if (FALSE == status)
  {
    error_printf("Error checking reactants.\n");
    return FALSE;
  }

#ifdef RULE_CHECK_DEBUG
  DYNPS_SPECIES_OUTPUT(OUTPUT_DEBUG, 4, dynps_species,
                       world->dynspecies, world->dyncolors);
#endif

  *pps_particles =
    dynps_species_setup_particles(world->dynspecies, pparticle_count,
                                  dynps_species, world->dyncolors);
  if (NULL == *pps_particles)
  {
    error_printf("Error setting up particles.\n");
    return FALSE;
  }

  return TRUE;
}

/* \brief Output all the current state information.
**
** \param world Current state
** \param log_time Current simulation time
*/
void
particlestoc_output_all(const output_type ptype, const int level,
                        world_t *world, const double log_time)
{
  // Print current time.
  custom_printf(ptype, level, "Current time: %e\n", log_time);

  dynspecies_output(ptype, level, world->dynspecies, world->dyncolors);
  reaction_definition_list_output(ptype, level,
                                  world->reaction_definition_list,
                                  world->dynrules, world->dyncolors);
  report_definition_list_output(ptype, level, world);
}

/* \brief Run the simulation as defined
**
** Precondition: World has been defined.
** Postcondition: Max_time has been reached, or an error occured.
**
** \param world State information
** \param logging_interval Amount of simulated time between log entries
** \param report_callback Callback function used to output reports
** \param max_time When to terminate simulation
** \param volume Volume simulation occurs in
** \param time_delta time step used for STOCHSIM algorithm
*/
boolean
particlestoc_run(world_t *world, const double orig_logging_interval,
                 report_callback_t report_callback,
                 const double max_time, const double volume,
                 const double orig_time_delta, const double max_prob,
                 const double pseudo_factor,
                 const double pseudo_count, const double avogadro)
{
  double time_delta;
  double logging_interval;
  double particle1_elements;
  double particle2_elements;
  boolean status;
  double next_time_to_log;
  size_t particle_count;
  double pseudo_particle_count;
  ps_particle_t *ps_particles;
  dynps_species_t *dynps_species;
  size_t particle1_index;
  size_t particle2_index;
  dynps_reaction_containers_t *dynps_reaction_containers;

  // Always log the starting condition.
  next_time_to_log = 0;

  dynps_reaction_containers =
    dynps_reaction_containers_create(INITIAL_REACTION_CONTAINERS);
  if (NULL == dynps_reaction_containers)
  {
    assert(0);
    return FALSE;
  }

  /* Sanity check */
  nauty_check(WORDSIZE, 1, 1, NAUTYVERSIONID);

  time_delta = orig_time_delta;

  if (orig_logging_interval < 0.0)
  {
    logging_interval = 0;
  }
  else
  {
    logging_interval = orig_logging_interval;
  }

  world->current_time = 0;
  ps_particles = NULL;

  // Create structure to store ps_specie information.
  dynps_species = dynps_species_create(INITIAL_DYNPS_SPECIES);

  // Finish setting up simulation.
  status = particlestoc_complete(world, dynps_species, &ps_particles,
                                 &particle_count);
  if (FALSE == status)
  {
    error_printf("Error setting up particlestoc environment.\n");
    return FALSE;
  }

  // Calculate the optimal number of pseudo-particles. */
  pseudo_particle_count =
    particlestoc_get_pseudo_particle_count(world->reaction_definition_list,
                                           particle_count, volume,
                                           pseudo_factor, pseudo_count,
                                           avogadro);
   if (0 == pseudo_particle_count)
   {
     error_printf("Error calculating number of pseudo-particles needed.\n");
     return FALSE;
   }

  // Update probabilities based on molecule counts and rate constants.
  // Will recompute time_delta if any probability greater than max_prob.
  ps_reaction_definition_init_all_prob(particle_count, pseudo_particle_count,
                                       &time_delta, volume,
                                       world->reaction_definition_list,
                                       max_prob, avogadro);

  // Want random number from 0 to (count - 1).
  particle1_elements = particle_count;
  particle2_elements = (double)particle_count +
                     (double)pseudo_particle_count;

  DEBUG_PRINTF(1, "Particles: %zu, Pseudo-particles: %.0f\n", particle_count,
               pseudo_particle_count);

  DEBUG_PRINTF(6, "p1_factor: %f, p2_factor: %f\n", particle1_elements,
               particle2_elements);

  DEBUG_PRINTF(6, "Species:\n");
  DYNPS_SPECIES_OUTPUT(OUTPUT_DEBUG, 6, dynps_species,
                       world->dynspecies, world->dyncolors);

  DYNCOLORS_OUTPUT(OUTPUT_DEBUG, 6,
                   world->dyncolors, "", "", TRUE, FALSE);
  DEBUG_PRINTF(6, "Reaction definitions:\n");
  REACTION_DEFINITION_LIST_OUTPUT(OUTPUT_DEBUG, 6,
                                  world->reaction_definition_list,
                                  world->dynrules, world->dyncolors);
  DEBUG_PRINTF(6, "Report definitions:\n");
  REPORT_DEFINITION_LIST_OUTPUT(OUTPUT_DEBUG, 6, world);

  DEBUG_PRINTF(9, "Particles:\n");
  PS_PARTICLE_OUTPUT_ALL(OUTPUT_DEBUG, 9, ps_particles,
                         particle_count, dynps_species);

  // Output column headers
  report_definition_output_names(OUTPUT_STANDARD, 0, world);

  for (;;)
  {
#if RULE_CHECK_DEBUG
    DEBUG_PRINTF (6, "***** Starting generation %.0f (time:%e).\n",
                  world->current_gen, world->current_time);
#endif

    if (world->current_time >= max_time)
    {
      // Output last time course and end simulation.
      report_callback(OUTPUT_STANDARD, 0, world, max_time);
      break;
    }

    // Check if current state should be logged.
    while ((world->current_time >= next_time_to_log) &&
           (0.0 != logging_interval))
    {
      // Report current status.
      if (NULL != report_callback)
      {
        report_callback(OUTPUT_STANDARD, 0, world, next_time_to_log);
      }
      else
      {
        //world_output_all(OUTPUT_STANDARD, 0, world, next_time_to_log);
      }

      next_time_to_log = next_time_to_log + logging_interval;

      if (next_time_to_log >= max_time)
      {
        break;
      }
    }


    // Choose first particle.
    particle1_index = MYRAND_INDEX(particle1_elements);

    // Choose second particle or pseudo-particle.
    particle2_index = MYRAND_INDEX(particle2_elements);

#if RULE_CHECK_DEBUG
    DEBUG_PRINTF(6, "particle1: %zu, particle2: %zu.\n", particle1_index,
                 particle2_index);
#endif

    // Try to perform a reaction with the particle(s).
    if (particle2_index < particle_count)
    {
      // Selected a normal particle, so attempt bimolecular reaction.
      if (particle1_index != particle2_index)
      {
#if RULE_CHECK_DEBUG
        DEBUG_PRINTF(12, "Particles:\n");
        PS_PARTICLE_OUTPUT_ALL(OUTPUT_DEBUG, 12, ps_particles,
                               particle_count, dynps_species);
#endif

        // Picked two unique particles.
        status = ps_specie_attempt_reaction(ps_particles + particle1_index,
                                            ps_particles + particle2_index,
                                            ps_particles, world, dynps_species,
                                            dynps_reaction_containers);
        if (FALSE == status)
        {
          error_printf("Error attempting bimolecular reaction");
          return FALSE;
        }
      }
    }
    else
    {
      // Selected a pseudo-particle, so attempt unimolecular reaction.
      status = ps_specie_attempt_reaction(ps_particles + particle1_index,
                                          NULL, ps_particles,
                                          world, dynps_species,
                                          dynps_reaction_containers);
      if (FALSE == status)
      {
        error_printf("Error attempting unimolecular reaction");
        return FALSE;
      }
    }

    // Reset reaction containers.
    dynps_reaction_containers_reset(dynps_reaction_containers);

#if RULE_CHECK_DEBUG
    DEBUG_PRINTF(12, "Particles:\n");
    PS_PARTICLE_OUTPUT_ALL(OUTPUT_DEBUG, 12,
                           ps_particles, particle_count, dynps_species);
#endif

    // Increment time.
    world->current_time += time_delta;
    world->current_gen++;
  }

  DEBUG_PRINTF(6, "Species:\n");
  DYNPS_SPECIES_OUTPUT(OUTPUT_DEBUG, 6, dynps_species,
                       world->dynspecies, world->dyncolors);
  DEBUG_PRINTF(6, "Report definitions:\n");
  REPORT_DEFINITION_LIST_OUTPUT(OUTPUT_DEBUG, 6, world);
#if RULE_CHECK_DEBUG
  DEBUG_PRINTF(8, "Particles:\n");
  PS_PARTICLE_OUTPUT_ALL(OUTPUT_DEBUG, 8, ps_particles,
                         particle_count, dynps_species);
#endif
  DEBUG_PRINTF(0, "Number of time steps: %.0f\n", world->current_gen - 1);
  DEBUG_PRINTF(0, "Number of completed reactions: %.0f\n",
               world->reaction_count);
  DEBUG_PRINTF(0, "Number of species: %zu\n", world->dynspecies->count);
  if (CUSTOM_OUTPUT_CHECK(5))
  {
    ilist_t *current_ilist;
    size_t unused_count;
    int max_size;
    specie_t *current_specie;
    int current_size;
    size_t i;

    // Count the currently unused species.
    unused_count = 0;
    current_ilist = dynps_species->available.a;
    while (current_ilist)
    {
      unused_count += 1;
      current_ilist = current_ilist->next;
    }
    DEBUG_PRINTF(0, "Unused species: %zu\n", unused_count);

    // Find the size of the largest graph.
    max_size = 0;
    for (i = 0, current_specie = world->dynspecies->species.a;
         i < world->dynspecies->count;
         i++, current_specie++)
    {
      current_size = current_specie->cg->n;

      if (current_size > max_size)
      {
        max_size = current_size;
      }
    }

    DEBUG_PRINTF(0, "Max graph size: %d\n", max_size);
  }

  return TRUE;
}


/* \brief Calculate pseudo-particle count
**
** \note Currently just sets pseudo_count equal to real_count.
**
** \param reaction_definition_list List of all defined reactions
** \param particle_count Number of particles in system
*/
static double
particlestoc_get_pseudo_particle_count
  (reaction_definition_t *reaction_definition_list,
   const size_t particle_count, const double volume,
   const double pseudo_factor,
   const double orig_pseudo_count, const double avogadro)
{
  double pseudo_count;

  if (orig_pseudo_count != 0)
  {
    pseudo_count = orig_pseudo_count;
  }
  else if (pseudo_factor != 0)
  {
    pseudo_count = particle_count * pseudo_factor;
  }
  else
  {
    double max_uni;
    double max_bi;


    // Use STOCHSIM algorithm for choosing pseudo particle count.

    reaction_definition_find_max_rc(reaction_definition_list,
                                    &max_uni, &max_bi);

    if ((0 == max_bi) && (0 == max_uni))
    {
      // Error, this shouldn't happen.
      error_printf("Max uni and bimolecular rate constant is zero.\n");
      assert(0);
      return 0;
    }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(2, "max_uni:%e, max_bi:%e\n", max_uni, max_bi);
#endif

    if (0 == max_uni)
    {
      // No unimolecular reactions, so pseudo particle not needed.
      pseudo_count = 1;
    }
    else if (0 == max_bi)
    {
      // No bimolecular reactions, so simply use large number of
      // pseudo particles.
      pseudo_count = particle_count * 1.0e5;
    }
    else
    {
#if 1
      pseudo_count = 2 * (max_uni/max_bi);

      // Round to the non-zero positive integer nearest to x.
      pseudo_count = floor(pseudo_count + 0.5);
#else
      pseudo_count = particle_count;
#endif

      // Should at least be 1.
      if (pseudo_count < 1)
      {
        pseudo_count = 1;
      }
    }
  }

  return pseudo_count;
}

