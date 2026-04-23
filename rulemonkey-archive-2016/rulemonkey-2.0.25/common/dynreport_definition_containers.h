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

#ifndef DYNREPORT_DEFINITION_CONTAINERS_H
#define DYNREPORT_DEFINITION_CONTAINERS_H

typedef struct dynreport_definition_containers_struct dynreport_definition_containers_t;

typedef union report_definition_containerarray_union
        report_definition_containerarray_t;

#include "dynrules.h"
#include "rule.h"
#include "report_definition_container.h"
#include "dyncolors.h"
#include "world.h"

// To prevent type-punning warnings.
union report_definition_containerarray_union
{
  report_definition_container_t *a;
  void *v;
};

struct dynreport_definition_containers_struct
{
  // Pointer to allocated array of report_definition_containers.
  report_definition_containerarray_t report_definition_containers;

  // Number of elements currently in use.
  size_t count;

  // Number of elements currently allocated.
  size_t alloc_count;
};

boolean dynreport_definition_containers_report_definition_container_create
  (report_definition_t *report_definition, const int pop,
   const int multiplier,
   dynreport_definition_containers_t *dynreport_definition_containers);
void dynreport_definition_containers_remove_count
       (dynreport_definition_containers_t *dynreport_definition_containers, world_t *world);
void dynreport_definition_containers_update_total
  (dynreport_definition_containers_t *dynreport_definition_containers,
   const int delta_total, world_t *world);
void dynreport_definition_containers_output
  (const output_type ptype, const int level,
   dynreport_definition_containers_t *dynreport_definition_containers);
report_definition_container_t *dynreport_definition_containers_report_definition_container_alloc
                      (dynreport_definition_containers_t *dynreport_definition_containers);
dynreport_definition_containers_t *dynreport_definition_containers_create(const size_t init_count);
void dynreport_definition_containers_destroy
       (dynreport_definition_containers_t *dynreport_definition_containers);

#endif
