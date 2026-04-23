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


#include "ps_reaction_definition.h"
#include "constants.h"
#include <assert.h>
#include "output.h"

void suppress_warnings33() {SUPPRESS_WARNINGS;}

/**
 * Computes probability that a given reaction will occur given that the
 * appropriate particles were randomly selected.  Will recompute time_delta
 * if any of the probabilities are greater than max_prob.
 */
void
ps_reaction_definition_init_all_prob
  (const double particle_count, const double pseudo_particle_count,
   double *ptime_delta, const double volume,
   reaction_definition_t *reaction_definition_list, const double max_prob,
   const double avogadro)
{
  reaction_definition_t *current_reaction_definition;
  double numerator;
  double factor1;
  double factor2;
  double largest_prob;

  largest_prob = 0;

  if (0 == *ptime_delta)
  {
    // Time delta not initialized, set to 1 initially.
    *ptime_delta = 1;
  }

  numerator =
    (double)particle_count *
     ((double)particle_count + (double)pseudo_particle_count) *
     *ptime_delta;
  if (0 == pseudo_particle_count)
  {
    factor1 = 0;
  }
  else
  {
    factor1 = numerator / (double)pseudo_particle_count;
  }

  factor2 = numerator / (double)2.0;

  current_reaction_definition = reaction_definition_list;
  while (NULL != current_reaction_definition)
  {
    if ((1 == current_reaction_definition->reactant_definition_count) &&
        (0 == current_reaction_definition->self_bind_count) &&
        (0 == current_reaction_definition->transient_self_bind_count))
    {
      // Unimolecular and uniparticle reaction.
      assert(factor1 > 0);
      current_reaction_definition->prob =
        current_reaction_definition->rate_const * factor1;
    }
    else
    {
      int self_bind_test;
      self_bind_test = current_reaction_definition->self_bind_count +
                       current_reaction_definition->transient_self_bind_count;

      if (((current_reaction_definition->reactant_definition_count > 2) ||
           (self_bind_test > 1)) ||
          ((current_reaction_definition->reactant_definition_count > 1) &&
           (self_bind_test != 0)))
      {
        error_printf("Do not support more than two reacting particles.\n");
        assert(0);
        exit(1);
      }
      // Bimolecular or biparticle reaction.
      current_reaction_definition->prob =
        current_reaction_definition->rate_const * factor2;
    }

    assert(current_reaction_definition->prob >= 0);

    if (current_reaction_definition->prob > largest_prob)
    {
      // Save largest probability in case we need to recompute time_delta.
      largest_prob = current_reaction_definition->prob;
    }

    DEBUG_PRINTF(2, "Reaction definition:%s, rc:%e, prob:%e\n",
                 current_reaction_definition->name,
                 current_reaction_definition->rate_const,
                 current_reaction_definition->prob);

    current_reaction_definition = current_reaction_definition->next;
  }

  if (largest_prob > max_prob)
  {
    double new_factor;

    // Get factor to update time delta and all probabilities.
    new_factor = max_prob / largest_prob;

    // Update time_delta.
    *ptime_delta *= new_factor;
    DEBUG_PRINTF(1, "New time delta:%e\n", *ptime_delta);

    // Update each probability.
    current_reaction_definition = reaction_definition_list;
    while (NULL != current_reaction_definition)
    {
      current_reaction_definition->prob *= new_factor;
      DEBUG_PRINTF(1, "Updated reaction definition:%s, rc:%e, prob:%e\n",
                   current_reaction_definition->name,
                   current_reaction_definition->rate_const,
                   current_reaction_definition->prob);

      current_reaction_definition = current_reaction_definition->next;
    }
  }
  DEBUG_PRINTF(5, "Time delta:%e\n", *ptime_delta);
}
