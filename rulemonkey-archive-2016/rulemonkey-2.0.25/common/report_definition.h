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

#ifndef REPORT_DEFINITION_H
#define REPORT_DEFINITION_H

typedef struct report_definition_struct report_definition_t;

#include "llist.h"
#include "dynrules.h"
#include "output.h"

#define REPORT_DEFINITION_OUTPUT_ALL_REPORTS(ptype, level, world, log_time) \
  CHECKDEBUG(report_definition_output_all_reports, ptype, level, \
             world, log_time)

#define REPORT_DEFINITION_LIST_OUTPUT(ptype, level, world) \
  CHECKDEBUG(report_definition_list_output, ptype, level, world)

typedef void(*report_callback_t)(const output_type ptype, const int level,
             world_t *world, double time);

struct report_definition_struct
{
  // For reporting.
  const char *name;

  // Current count of all matched species.
  double total;

  // Rule to use.
  int rule_index;

  // Pointer to rule to use.
  rule_t *rule;

  // Number to divide matched with.
  double denominator;

  // Implement linked list of report definitions.
  report_definition_t *next;
};

boolean report_definition_check_all
          (specie_t *specie, dynrules_t *dynrules,
           report_definition_t *report_definition_list, world_t *world);
boolean report_definition_check(specie_t *specie, dynrules_t *dynrules,
                                report_definition_t *report_definition_list,
                                world_t *world);
boolean report_definition_add(const char *name, const int rule_index,
                              double denominator, world_t *world);
report_definition_t *report_definition_create
                       (const char *name, const int rule_index,
                        double denominator, report_definition_t *next);

void report_definition_output_names (const output_type ptype, const int level,
                                     world_t *world);
void report_definition_output_all_reports
  (const output_type ptype, const int level, world_t *world, double log_time);
void report_definition_list_output(const output_type ptype, const int level,
                                   world_t *world);
boolean report_definition_complete(report_definition_t *report_definition,
                                   dynrules_t *dynrules);

#endif
