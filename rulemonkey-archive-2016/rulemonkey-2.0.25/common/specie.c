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

#include "specie.h"
#include "naututil.h"
#include "constants.h"
#include "output.h"

void suppress_warnings10() {SUPPRESS_WARNINGS;}

boolean
specie_setup(int id, cgraph_t *cg, const char *name, const double pop,
             boolean update, boolean preserve, specie_t *new_specie)
{
  int node_index;
  int adjacent_count;
  set *row;
  int i;

  // Set structure members.
  new_specie->id = id;
  if (new_specie->name != NULL)
  {
    free(new_specie->name);
  }
  if (name != NULL)
  {
    new_specie->name = strdup(name);
  }
  else
  {
    new_specie->name = NULL;
  }
  new_specie->pop = pop;
  new_specie->update = update;
  new_specie->preserve = preserve;
  if (new_specie->cg != NULL)
  {
    int status;

    status = cgraph_ncopy(cg, new_specie->cg);
    if (FALSE == status)
    {
      free(new_specie->name);
      return FALSE;
    }
  }
  else
  {
    new_specie->cg = cgraph_copy(cg);
    if (NULL == new_specie->cg)
    {
      assert(0);
      free(new_specie);
      return FALSE;
    }
  }

  if (new_specie->dynreport_definition_containers != NULL)
  {
    new_specie->dynreport_definition_containers->count = 0;
  }
  else
  {
    new_specie->dynreport_definition_containers = 
      dynreport_definition_containers_create
        (INITIAL_DYNREPORT_DEFINITION_CONTAINERS);
  }

  // For each node, compute the number of adjacent nodes.
  for (node_index = 0; node_index < cg->n; node_index++)
  {
    adjacent_count = 0;

    row = GRAPHROW(cg->g.a, node_index, cg->m);

    for (i = -1; (i = nextelement(row, cg->m, i)) >= 0;)
    {
      adjacent_count++;
    }

    new_specie->cg->adjacent_counts.a[node_index] = adjacent_count;
  }

  return TRUE;
}

void
specie_update(specie_t *specie, const double pop_delta, world_t *world)
{
  if (TRUE == specie->update)
  {
    specie->pop += pop_delta;

    dynreport_definition_containers_update_total
     (specie->dynreport_definition_containers, pop_delta, world);
  }

#ifndef DS_NO_DOUBLECHECK
  assert(specie->pop >= 0);
#endif
}

int
specie_check(specie_t *specie, dynrule_aggregates_t *dynrule_aggregates,
             reaction_definition_t *reaction_definition, int reactant_index,
             dynrules_t *dynrules, dyncolors_t *dyncolors, int *pcount)
{
  reactant_definition_t *reactant_definition;
  int status;

  reactant_definition =
    &reaction_definition->reactant_definitions[reactant_index];

#ifdef RULE_CHECK_DEBUG
  DEBUG_PRINTF(8, "About to call rule_check.\n");
  DEBUG_PRINTF(8, "specie:.\n");
  SPECIE_OUTPUT(OUTPUT_DEBUG, 8, specie, dyncolors);
#endif
  status = rule_check(specie, dynrule_aggregates,
                      reactant_definition->rule_list,
                      pcount, dynrules,
                      dyncolors);
  if (0 == status)
  {
    error_printf("Error with reaction definition: %s\n",
                 reaction_definition->name);
    return 0;
  }

  // No error occured.
  return 1;
}

void
specie_output(const output_type ptype, const int level,
              specie_t *specie, dyncolors_t *dyncolors)
{
  CUSTOM_PRINTF(ptype, level, "Specie: %d", specie->id + 1);
  if (NULL != specie->name)
  {
    CUSTOM_PRINTF(ptype, level, " (\"%s\")", specie->name);
  }
  CUSTOM_PRINTF(ptype, level, "\n");
  CUSTOM_PRINTF(ptype, level, "Population: %.0f\n", specie->pop);
  CUSTOM_PRINTF(ptype, level, "Graph:\n");

  cgraph_output_bngl(ptype, level, specie->cg, dyncolors);
  CUSTOM_PRINTF(ptype, level, "\n");

  dynreport_definition_containers_output
    (ptype, level, specie->dynreport_definition_containers);
}

void
specie_output_net(const output_type ptype, const int level,
                  specie_t *specie, dyncolors_t *dyncolors)
{
  int i;
  color_t *color;

  CUSTOM_PRINTF(ptype, level, "%d (", specie->id);
  for (i = 0; i < specie->cg->n; i++)
  {
    if (i != 0)
    {
      CUSTOM_PRINTF(ptype, level, ",");
    }
    color = dyncolors_color_value_lookup(specie->cg->colors.a[i], dyncolors);
    if (NULL != color)
    {
      CUSTOM_PRINTF(ptype, level, "%s", color->name);
    }
    else
    {
      // Unable to find color name.
      CUSTOM_PRINTF(ptype, level, "%d", specie->cg->colors.a[i]);
    }
  }
  CUSTOM_PRINTF(ptype, level, ") %.0f", specie->pop);
  CUSTOM_PRINTF(ptype, level, "\n");
}

