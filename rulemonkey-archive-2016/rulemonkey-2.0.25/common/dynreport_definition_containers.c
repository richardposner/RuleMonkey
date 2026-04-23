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
#include "dynreport_definition_containers.h"
#include "dynarray.h"
#include "constants.h"
#include "output.h"

void suppress_warnings18() {SUPPRESS_WARNINGS;}

// Local function declarations.
static int dynreport_definition_containers_resize
  (dynreport_definition_containers_t *dynreport_definition_containers,
   const size_t count);

boolean
dynreport_definition_containers_report_definition_container_create
  (report_definition_t *report_definition, const int pop,
   const int multiplier,
   dynreport_definition_containers_t *dynreport_definition_containers)
{
  report_definition_container_t *new_report_definition_container;

  // Need to create new report_definition_container.
  new_report_definition_container =
    dynreport_definition_containers_report_definition_container_alloc
      (dynreport_definition_containers);
  if (NULL == new_report_definition_container)
  {
    assert(0);
    return FALSE;
  }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "UPDATE [NEW] (%s) %f += %d * %d\n",
                 report_definition->name,
                 report_definition->total,
                 multiplier, pop);
#endif

  // Add initial population.
  report_definition->total += multiplier * pop;

  // Report definitions are stored in linked list, so ptr is safe to use.
  new_report_definition_container->report_definition = report_definition;
  new_report_definition_container->multiplier = multiplier;

  return TRUE;
}

void
dynreport_definition_containers_update_total
  (dynreport_definition_containers_t *dynreport_definition_containers,
   const int delta_total, world_t *world)
{
  int i;
  int report_definition_containers_count;
  report_definition_container_t *current_report_definition_container;
  report_definition_t *current_report_definition;


  report_definition_containers_count = dynreport_definition_containers->count;
  for (i = 0, current_report_definition_container =
                dynreport_definition_containers->report_definition_containers.a;
       i < report_definition_containers_count;
       i++, current_report_definition_container++)
  {
    current_report_definition =
      current_report_definition_container->report_definition;

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "UPDATE (%s) %f += %d * %d\n",
                 current_report_definition->name,
                 current_report_definition->total,
                 current_report_definition_container->multiplier,
                 delta_total);
#endif

    current_report_definition->total +=
      current_report_definition_container->multiplier * delta_total;
  }
}

void
dynreport_definition_containers_output
  (const output_type ptype, const int level,
   dynreport_definition_containers_t *dynreport_definition_containers)
{
  int i;
  int report_definition_containers_count;
  report_definition_container_t *current_report_definition_container;

  CUSTOM_PRINTF(ptype, level, "Report Definition Containers:\n");

  report_definition_containers_count = dynreport_definition_containers->count;
  if (0 == report_definition_containers_count)
  {
    custom_printf(ptype, level, "  None\n");
    return;
  }
  for (i = 0, current_report_definition_container =
                dynreport_definition_containers->report_definition_containers.a;
       i < report_definition_containers_count;
       i++, current_report_definition_container++)
  {
    CUSTOM_PRINTF
      (ptype, level, "  pop * %d : %s\n",
       current_report_definition_container->multiplier,
       current_report_definition_container->report_definition->name);
  }
}
   
   

report_definition_container_t *
dynreport_definition_containers_report_definition_container_alloc
  (dynreport_definition_containers_t *dynreport_definition_containers)
{
  int status;

  status =
    dynreport_definition_containers_resize
      (dynreport_definition_containers,
       dynreport_definition_containers->count + 1);
  if (FALSE == status)
  {
    assert(0);
    return NULL;
  }

  return &dynreport_definition_containers->report_definition_containers.a[dynreport_definition_containers->count++];
}

dynreport_definition_containers_t *
dynreport_definition_containers_create(const size_t init_count)
{
  dynreport_definition_containers_t *new_dynreport_definition_containers;
  int status;

  // All fields need to be initialized to zero anyways, so use calloc.
  new_dynreport_definition_containers = calloc(1, sizeof(dynreport_definition_containers_t));
  if (NULL == new_dynreport_definition_containers)
  {
    error_printf("Unable to allocate new dynreport_definition_containers.\n");
    return NULL;
  }

  status =
  dynarray_create
    (init_count,
     &new_dynreport_definition_containers->alloc_count,
     sizeof(report_definition_container_t),
     &new_dynreport_definition_containers->
              report_definition_containers.v);
  if (0 == status)
  {
    // Error message output in dynarray_create.
    free(new_dynreport_definition_containers);
    return NULL;
  }

  return new_dynreport_definition_containers;
}

static int
dynreport_definition_containers_resize
  (dynreport_definition_containers_t *dynreport_definition_containers,
   const size_t count)
{
  int status;

  status =
    dynarray_resize
      (count,
       &dynreport_definition_containers->alloc_count,
       sizeof(report_definition_container_t),
       &dynreport_definition_containers->report_definition_containers.v);
  if (0 == status)
  {
    // Error message output in dynarray_resize.
    return 0;
  }

  return 1;
}

void
dynreport_definition_containers_destroy
  (dynreport_definition_containers_t *dynreport_definition_containers)
{
  int i;
  report_definition_container_t *current_report_definition_container;

  for (i = 0,
        current_report_definition_container = dynreport_definition_containers->report_definition_containers.a;
       i <dynreport_definition_containers->count;
       i++, current_report_definition_container++)
  {
    //report_definition_container_destroy(current_report_definition_container);
  }

  free(dynreport_definition_containers->report_definition_containers.a);
  free(dynreport_definition_containers);
}
