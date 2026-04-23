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
#include "world.h"
#include "constants.h"
#include "output.h"

void suppress_warnings14() {SUPPRESS_WARNINGS;}

world_t *
world_create()
{
  world_t *new_world;
  double d;
  double a;

  // Most attributes need to be initialized to 0 anyways.
  new_world = calloc(1, sizeof(world_t));
  if (NULL == new_world)
  {
    perror("Allocating new world");
    return NULL;
  }

  new_world->dynrules = dynrules_create(INITIAL_DYNRULES);
  if (NULL == new_world->dynrules)
  {
    perror("Allocating dynrules");
    return NULL;
  }
  new_world->dyncolors = dyncolors_create(INITIAL_DYNCOLORS);
  if (NULL == new_world->dyncolors)
  {
    perror("Allocating dyncolors");
    return NULL;
  }
  new_world->dynspecies = dynspecies_create(INITIAL_DYNSPECIES);
  if (NULL == new_world->dynspecies)
  {
    perror("Allocating dynspecies");
    return NULL;
  }

  // Used when a reaction has REACTION_FLAG_DIFFUSION set.
  // Diffusion coeff: 2*10^-7.
  d = 2e-7;
  // 4*10^-4
  a = 4e-4;
  new_world->diffusion_denominator = 4 * 3.14159265 * d * a;

  return new_world;
}

boolean
world_complete(world_t *world)
{
  int status;
  reaction_definition_t *current_reaction_definition;
  report_definition_t *current_report_definition;

  // Resolve rule names to pointers.
  status = dynrules_complete(world->dynrules);
  if (0 == status)
  {
    // Error message already output by dynrules_complete.
    assert(0);
    return FALSE;
  }

  // Finish setting up report definitions.
  current_report_definition = world->report_definition_list;
  while (NULL != current_report_definition)
  {
    status = report_definition_complete(current_report_definition,
                                        world->dynrules);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }

    current_report_definition = current_report_definition->next;
  }

  // Finish setting up reaction definitions.
  current_reaction_definition = world->reaction_definition_list;
  while (NULL != current_reaction_definition)
  {
    status = reaction_definition_complete(current_reaction_definition,
                                          world->dynrules);
    if (0 == status)
    {
      assert(0);
      return FALSE;
    }

    current_reaction_definition = current_reaction_definition->next;
  }

  return TRUE;
}

void
world_output_pop(const output_type ptype, const int level,
                 world_t *world, double log_time)
{
  // Print current time.
  CUSTOM_PRINTF(ptype, level, "%e", world->current_time);

  dynspecies_output_pop(ptype, level, world->dynspecies, world);
}

