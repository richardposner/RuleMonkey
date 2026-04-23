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
#include <float.h>
#include "constants.h"
#include "reactionstoc.h"
#include "rs_reaction_definition.h"
#include "rs_specie_container.h"
#include "fsa.h"
#include "world.h"
#include "myrand.h"

void suppress_warnings41() {SUPPRESS_WARNINGS;}

static int reactionstoc_complete
             (fsa_t *rs_specie_fsa,
              rs_reaction_definition_t *rs_reaction_definitions,
              size_t reaction_definition_count,
              world_t *world, double *ptotal_prob);

//RG: Added by Ryan Gutenkunst (RG) to support continued simulations.
//    It might be better to store this within the world instance, but using a
//    global variable here minimizes changes that need to be made to the code.
int num_species_in_last_run = 0;

/* \brief Run the simulation as defined, use network-free
**  Gillespie direct method
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
int
reactionstoc_run(world_t *world, const double orig_logging_interval,
                 report_callback_t report_callback, const double max_time,
                 const double pseudo_factor, const double pseudo_count)
{
  boolean status;
  double logging_interval;
  double next_time_to_log;
  double time_delta;
  double total_prob;
  size_t reaction_definition_count;

  rs_reaction_definition_t *rs_reaction_definitions;
  fsa_t *rs_specie_fsa;

  // Always log the starting condition.
  //RG: But for continued simulations the starting condition is not necessarily
  //    at t = 0.
  next_time_to_log = world->current_time;

  /* Sanity check */
  nauty_check(WORDSIZE, 1, 1, NAUTYVERSIONID);

 if (orig_logging_interval < 0.0)
  {
    logging_interval = 0;
  }
  else
  {
    logging_interval = orig_logging_interval;
  }

  //RG: For continued simulations, we don't want to reset world time.
  //world->current_time = 0;

  // Count the number of reaction definitions.
  reaction_definition_count = 0;
  {
    reaction_definition_t *current_reaction_definition;

    current_reaction_definition = world->reaction_definition_list;
    while(current_reaction_definition != NULL)
    {
      reaction_definition_count++;
      current_reaction_definition = current_reaction_definition->next;
    }
  }

  // Create rs_specie_fsa.
  rs_specie_fsa = RS_SPECIES_CREATE();

  // Create rs_reaction_definitions.
  rs_reaction_definitions =
    calloc(reaction_definition_count, sizeof(rs_reaction_definition_t));

  // Complete any remaining computations needed before starting simulation.
  total_prob = 0;
  status = reactionstoc_complete(rs_specie_fsa, rs_reaction_definitions,
                                 reaction_definition_count,
                                 world, &total_prob);
  if (0 == status)
  {
    error_printf("Error setting up reactionstoc\n");
    return 0;
  }

  DEBUG_PRINTF(6, "Species:\n");
  RS_SPECIE_FSA_OUTPUT(OUTPUT_DEBUG, 6, rs_specie_fsa, world);

  DYNCOLORS_OUTPUT(OUTPUT_DEBUG, 6,
                   world->dyncolors, "", "", TRUE, FALSE);
  DEBUG_PRINTF(6, "Reaction definitions:\n");
  if (CUSTOM_OUTPUT_CHECK(6))
  {
    int i;
    rs_reaction_definition_t *current_rs_reaction_definition;

    for (i = 0, current_rs_reaction_definition = rs_reaction_definitions;
         i < reaction_definition_count;
         i++, current_rs_reaction_definition++)
    {
      RS_REACTION_DEFINITION_OUTPUT(OUTPUT_DEBUG, 6,
                                    current_rs_reaction_definition,
                                    world->dynrules, world->dyncolors);
    }
  }
  DEBUG_PRINTF(6, "Report definitions:\n");
  REPORT_DEFINITION_LIST_OUTPUT(OUTPUT_DEBUG, 6, world);

  //RG: Only output column headers if this is the first call to simulate_rs
  if(world->current_time == 0){
      // Output column headers
      report_definition_output_names(OUTPUT_STANDARD, 0, world);

      bng_standard_printf("# ");
      report_definition_output_names(OUTPUT_BNG_STANDARD, 0, world);
  }

  // Keep performing reactions until simulation ends.
  for (;;)
  {
    double current_time;

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF (6, "***** Starting generation %.0f (time:%e).\n",
                  world->current_gen, world->current_time);
    DEBUG_PRINTF(15, "Species:\n");
    RS_SPECIE_FSA_OUTPUT(OUTPUT_DEBUG, 6, rs_specie_fsa, world);
    DEBUG_PRINTF(6, "Reaction definitions:\n");
    if (CUSTOM_OUTPUT_CHECK(6))
    {
      int i;
      rs_reaction_definition_t *current_rs_reaction_definition;

      for (i = 0, current_rs_reaction_definition = rs_reaction_definitions;
           i < reaction_definition_count;
           i++, current_rs_reaction_definition++)
      {
        RS_REACTION_DEFINITION_OUTPUT(OUTPUT_DEBUG, 6,
                                      current_rs_reaction_definition,
                                      world->dynrules, world->dyncolors);
      }
    }
    else if (CUSTOM_OUTPUT_CHECK(2))
    {
      int i;
      rs_reaction_definition_t *current_rs_reaction_definition;

      custom_printf(OUTPUT_DEBUG, 2, "\n");
      for (i = 0, current_rs_reaction_definition = rs_reaction_definitions;
           i < reaction_definition_count;
           i++, current_rs_reaction_definition++)
      {
        rs_reaction_definition_short_output(OUTPUT_DEBUG, 2,
                                      current_rs_reaction_definition,
                                      world->dynrules, world->dyncolors);
      }
    }
    DEBUG_PRINTF(15, "Report definitions:\n");
    REPORT_DEFINITION_LIST_OUTPUT(OUTPUT_DEBUG, 6, world);
#endif

    /* It seems simpler to just set time_delta to be very,very large in this
     * case, so the appropriate logging still happens.
    if (0 == total_prob)
    {
      // Nothing to do.
      report_callback(OUTPUT_STANDARD, 0, world, max_time);

      // BioNetGen compatible output.
      bng_standard_printf(" ");
      report_callback(OUTPUT_BNG_STANDARD, 0, world, max_time);
      error_printf("Total probability is 0, nothing to do\n");
      break;
    }
    */

    //RG: Moved calculation of next reaction time to the front of the loop to
    //    ensure that state of system gets logged *before* the next reaction
    //    fires.
    // Get time to next reaction.
    if(total_prob > 0){
        time_delta = (1 / total_prob) * log(1 / MYRAND());
    }
    else{
        time_delta = DBL_MAX;
    }

    //RG: Set current time to be time of next reaction, so all the logging
    //    code works correctly.
    current_time = world->current_time + time_delta;

    // Check if current state should be logged.
    while ((current_time >= next_time_to_log) &&
           (0.0 != logging_interval))
    {
      // Report current status.
      if (NULL != report_callback)
      {
        report_callback(OUTPUT_STANDARD, 0, world, next_time_to_log);

        // BioNetGen compatible output.
        bng_standard_printf(" ");
        report_callback(OUTPUT_BNG_STANDARD, 0, world, next_time_to_log);
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

    if (current_time >= max_time)
    {
      // Output last time course and end simulation.
      report_callback(OUTPUT_STANDARD, 0, world, max_time);

      // BioNetGen compatible output.
      bng_standard_printf(" ");
      report_callback(OUTPUT_BNG_STANDARD, 0, world, max_time);
      //RG: We're breaking out of reaction loop, record the time we ended at.
      world->current_time = max_time;
      break;
    }

    //RG: Note that reaction will not be fired if it's firing time is greater
    //    than the maximum time requested in the simulation.
    status = rs_reaction_definition_perform_reaction
      (rs_reaction_definitions, reaction_definition_count,
       rs_specie_fsa, world, &total_prob);
    if (0 == status)
    {
      error_printf("Error attempting reaction\n");
      assert(0);
    }
    world->reaction_count++;

    // Increment time.
    world->current_time += time_delta;
    world->current_gen++;
  }

  DEBUG_PRINTF(5, "Species:\n");
  RS_SPECIE_FSA_OUTPUT(OUTPUT_DEBUG, 5, rs_specie_fsa, world);
  DEBUG_PRINTF(5, "Reaction definitions:\n");
  if (CUSTOM_OUTPUT_CHECK(5))
  {
    int i;
    rs_reaction_definition_t *current_rs_reaction_definition;

    for (i = 0, current_rs_reaction_definition = rs_reaction_definitions;
         i < reaction_definition_count;
         i++, current_rs_reaction_definition++)
    {
      RS_REACTION_DEFINITION_OUTPUT(OUTPUT_DEBUG, 5,
                                    current_rs_reaction_definition,
                                    world->dynrules, world->dyncolors);
    }
  }
  else if (CUSTOM_OUTPUT_CHECK(2))
  {
    int i;
    rs_reaction_definition_t *current_rs_reaction_definition;

    for (i = 0, current_rs_reaction_definition = rs_reaction_definitions;
         i < reaction_definition_count;
         i++, current_rs_reaction_definition++)
    {
      rs_reaction_definition_short_output(OUTPUT_DEBUG, 2,
                                    current_rs_reaction_definition,
                                    world->dynrules, world->dyncolors);
    }
  }
  DEBUG_PRINTF(5, "Report definitions:\n");
  REPORT_DEFINITION_LIST_OUTPUT(OUTPUT_DEBUG, 5, world);
  DEBUG_PRINTF(0, "Number of time steps: %.0f\n", world->current_gen);
  DEBUG_PRINTF(0, "Number of completed reactions: %.0f\n",
               world->reaction_count);
  DEBUG_PRINTF(0, "Number of species: %zu\n", world->dynspecies->count);
  if (CUSTOM_OUTPUT_CHECK(5))
  {
    llist_t *current_llist;
    size_t unused_count;
    int max_size;
    specie_t *current_specie;
    int current_size;
    size_t i;

    // Count the currently unused species.
    unused_count = 0;
    current_llist = rs_specie_fsa->available;
    while (current_llist != NULL)
    {
      unused_count += 1;
      current_llist = current_llist->next;
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

  //RG: Record total species in last run to support restarting simulations
  num_species_in_last_run = world->dynspecies->count;

  return TRUE;
}

static int
reactionstoc_complete(fsa_t *rs_specie_fsa,
                      rs_reaction_definition_t *rs_reaction_definitions,
                      size_t reaction_definition_count,
                      world_t *world, double *ptotal_prob)
{
  int status;
  reaction_definition_t *current_reaction_definition;
  rs_reaction_definition_t *current_rs_reaction_definition;
  int species_count;
  int i;
  specie_t *specie;

  status = world_complete(world);
  if (FALSE == status)
  {
    error_printf("Error setting up world\n");
    return FALSE;
  }

  current_reaction_definition = world->reaction_definition_list;
  current_rs_reaction_definition = rs_reaction_definitions;
  while(current_reaction_definition != NULL)
  {
    status = rs_reaction_definition_init(current_rs_reaction_definition,
                                         current_reaction_definition,
                                         world->dyncolors, world->dynrules);
    if (0 == status)
    {
      error_printf("Error running rs_reaction_definition_init\n");
      return 0;
    }

    current_rs_reaction_definition++;
    current_reaction_definition = current_reaction_definition->next;
  }

  // For each specie, find eligible reactant definitions.
  species_count = world->dynspecies->count;
  for (i = 0, specie = world->dynspecies->species.a;
       i < species_count;
       i++, specie++)
  {
    int status;
    rs_specie_t *rs_specie;


    // Create empty rs_specie for every specie.
    rs_specie = fsa_get_new(rs_specie_fsa);

    rs_specie->specie_id = i;

    //RG: Only update reports (observables) for species that have been added
    //    since the last run. Otherwise, we end up double counting all our
    //    previously observed species.
    if(i >= num_species_in_last_run){
        // Keep a note of which reports new specie should be included in.
        status = report_definition_check_all
            (specie, world->dynrules, world->report_definition_list, world);
        if (FALSE == status)
        {
          error_printf("Error running report_definition_check_all\n");
          return FALSE;
        }
    }

    // Find reactions that current specie matches.
    status = rs_specie_check_all_reactants(rs_specie, rs_reaction_definitions,
                                           specie, world, ptotal_prob);
    if (0 == status)
    {
      error_printf("Error running rs_specie_check_all_reactants\n");
      return 0;
    }
  }

  return 1;
}

