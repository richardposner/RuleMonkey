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

#include "report_definition.h"
#include "constants.h"
#include "output.h"
#include "rule.h"

void suppress_warnings16() {SUPPRESS_WARNINGS;}

/**
 * @brief Loop through every report definition to find all matches.
 *
 * @param specie Specie to check against every report definition.
 * @param dynrules All defined rules.
 * @param report_definition_list All report definitions.
 * @param world All data on current simulation state.
 * @returns FALSE if failure, TRUE otherwise.
 */
boolean
report_definition_check_all(specie_t *specie, dynrules_t *dynrules,
                            report_definition_t *report_definition_list,
                            world_t *world)
{
  report_definition_t *current_report_definition;
  boolean status;

  // Loop through every report_definition to find all matches.
  for (current_report_definition = report_definition_list;
       NULL != current_report_definition;
       current_report_definition = current_report_definition->next)
  {
    status = report_definition_check(specie, dynrules,
                                     current_report_definition, world);
    if (FALSE == status)
    {
      // Error occured.
      assert(0);
      return FALSE;
    }
  }

  return TRUE;
}

/**
 * @brief Check individual report definition against given specie.
 *
 * @param specie Specie to check against.
 * @param dynrules All defined rules.
 * @param report_definition Report definition to check specie against.
 * @param world All data on current simulation state.
 * @returns FALSE if failure, TRUE otherwise.
 */
boolean
report_definition_check(specie_t *specie, dynrules_t *dynrules,
                        report_definition_t *report_definition, world_t *world)
{
  int status;
  int adjusted_count;
  int count;

  status = rule_check(specie, NULL, report_definition->rule,
                      &count, dynrules, world->dyncolors);
  if (0 == status)
  {
    error_printf("Error with report definition: %s", report_definition->name);
    return FALSE;
  }

  if (count > 0)
  {
    if (0 == report_definition->denominator)
    {
      // Denominator of 0 means ignore count.
      adjusted_count = 1;
    }
    else
    {
      adjusted_count = count / report_definition->denominator;
    }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(3,
                 "Report def: %s, specie %d, count %d, adjusted count %d\n",
                 report_definition->name, specie->id, count,
                 adjusted_count);
#endif

    // Rule matched, add report definition to list.
    status =
      dynreport_definition_containers_report_definition_container_create
        (report_definition, specie->pop, adjusted_count,
         specie->dynreport_definition_containers);
    if (FALSE == status)
    {
      perror("report_definition_check failed to create new "
             "report definition container");
      assert(0);
      return FALSE;
    }
  }

  return TRUE;
}


/**
 * @brief Create definition and add it to the world list.
 *
 * @param name Name of new report.
 * @param rule_index Rule used to define report definition.
 * @param denominator Constant to divide result with.
 * @param world All data on current simulation state.
 */
boolean
report_definition_add(const char *name, const int rule_index,
                      double denominator, world_t *world)
{
  report_definition_t *new_report_definition;

  new_report_definition =
    report_definition_create(name, rule_index, denominator, NULL);
  if (NULL == new_report_definition)
  {
    assert(0);
    return FALSE;
  }

  if (NULL == world->report_definition_last)
  {
    // First report definition.
    world->report_definition_list = new_report_definition;
    world->report_definition_last = new_report_definition;
  }
  else
  {
    // Add to end of existing report definitions.
    world->report_definition_last->next = new_report_definition;
    world->report_definition_last = new_report_definition;
  }

  return TRUE;
}

/**
 * @brief Create definition
 *
 * @param name Name of new report.
 * @param rule_index Rule used to define report definition.
 * @param denominator Constant to divide result with.
 * @param world All data on current simulation state.
 * @return report definition object.
 */
report_definition_t *
report_definition_create(const char *name, const int rule_index,
                         double denominator, report_definition_t *next)
{
  report_definition_t *new_report_definition;

  new_report_definition = calloc(1, sizeof(report_definition_t));
  if (NULL == new_report_definition)
  {
    perror("Allocating new report definition");
    assert(0);
    return NULL;
  }

  // Set structure members.
  new_report_definition->name = strdup(name);
  if (NULL == new_report_definition->name)
  {
    perror("Allocating new report definition->name");
    free(new_report_definition);
    assert(0);
    return NULL;
  }

  new_report_definition->rule_index = rule_index;
  new_report_definition->denominator = denominator;
  new_report_definition->next = next;

  return new_report_definition;
}

/**
 * @brief Output the name of each report, used to define the column names.
 */
void
report_definition_output_names(const output_type ptype, const int level,
                               world_t *world)
{
  report_definition_t *current_report_definition;

  // Header for time.
  CUSTOM_PRINTF(ptype, level, "%-14s", "time");

  current_report_definition = world->report_definition_list;
  while (NULL != current_report_definition)
  {
    // Make sure there is always at least one space.
    CUSTOM_PRINTF(ptype, level, " %14s", current_report_definition->name);

    current_report_definition = current_report_definition->next;
  }

  CUSTOM_PRINTF(ptype, level, "\n");
}

/**
 * @brief Output the value for each report, used for column values.
 */
void
report_definition_output_all_reports(const output_type ptype, const int level,
                                     world_t *world, double log_time)
{
  report_definition_t *current_report_definition;

  CUSTOM_PRINTF(ptype, level, "%.8e", log_time);

  current_report_definition = world->report_definition_list;
  while (NULL != current_report_definition)
  {
    CUSTOM_PRINTF(ptype, level, " %.8e", current_report_definition->total);

    current_report_definition = current_report_definition->next;
  }

  CUSTOM_PRINTF(ptype, level, "\n");
}

/**
 * @brief Output all information about each report, used for debugging.
 */
void
report_definition_list_output(const output_type ptype, const int level,
                              world_t *world)
{
  report_definition_t *current_report_definition;

  current_report_definition = world->report_definition_list;
  while (NULL != current_report_definition)
  {
    custom_printf(ptype, level, "Observable: %s\n",
                  current_report_definition->name);
    custom_printf(ptype, level, "  Total: %.8e\n",
                  current_report_definition->total);
    custom_printf(ptype, level, "  Denominator: %f\n",
                  current_report_definition->denominator);

    dynrules_output_rule(ptype, level,
                         current_report_definition->rule, world->dynrules,
                         world->dyncolors);


    current_report_definition = current_report_definition->next;
  }

  custom_printf(ptype, level, "\n");
}

boolean
report_definition_complete(report_definition_t *report_definition,
                           dynrules_t *dynrules)
{
  report_definition->rule = dynrules->rules.a + report_definition->rule_index;

  return 1;
}

