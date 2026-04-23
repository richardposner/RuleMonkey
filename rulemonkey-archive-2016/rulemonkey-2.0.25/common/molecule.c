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
#include "molecule.h"
#include "dyncolors.h"

void suppress_warnings11() {SUPPRESS_WARNINGS;}

// Local function declarations.
static boolean molecule_component_resize(molecule_t *molecule,
                                         const size_t count);

boolean
molecule_setup(const int color_value, const size_t initial_component_count,
               molecule_t *new_molecule)
{
  int status;

  new_molecule->color_value = color_value;

  // Initialize counts.
  new_molecule->component_count = 0;
  new_molecule->component_alloc_count = 0;

  status = dynarray_create3(initial_component_count,
                            &new_molecule->component_alloc_count,
                            sizeof(int),
                            &new_molecule->component_color_values.v,
                            sizeof(int),
                            &new_molecule->component_multipliers.v,
                            sizeof(dynints_t *),
                            &new_molecule->state_dynints_array.v);
  if (0 == status)
  {
    // Error message output in dynarray_create3.
    return FALSE;
  }

  return TRUE;
}

int *
molecule_component_multipliers_table_create(molecule_t *molecule)
{
  return calloc(1, molecule->component_count * sizeof(int));
}

boolean
molecule_component_verify(const int component_color_value,
                          const int component_multiplier,
                          dynints_t *state_dynints,
                          int *molecule_component_multipliers_table,
                          dyncolors_t *dyncolors,
                          molecule_t *molecule,
                          int molecule_types_present)
{
  int component_count;
  int i;
  int j;
  int *current_component_color_value;

  // Find component index.
  component_count = molecule->component_count;
  for (i = 0,
       current_component_color_value = molecule->component_color_values.a;
       i < component_count;
       i++, current_component_color_value++)
  {
    if ((component_color_value == *current_component_color_value) ||
        (COLOR_WILDCARD == *current_component_color_value))
    {
      molecule_component_multipliers_table[i] += component_multiplier;
      if (molecule_component_multipliers_table[i] >
          molecule->component_multipliers.a[i])
      {
        // Too many components are present in the current molecule definition.
        return FALSE;
      }

      if (NULL != state_dynints)
      {
        for (j = 0; j < state_dynints->count; j++)
        {
          int state_color;

          state_color = state_dynints->ints.a[j];

          if ((COLOR_WILDCARD != state_color) &&
              (NULL != molecule->state_dynints_array.a[i]) &&
              (0 == dynints_int_exists(state_dynints->ints.a[j],
                                       molecule->state_dynints_array.a[i])))
          {
            if (1 == molecule_types_present)
            {
              error_printf
                ("State '%s' not defined for %s(%s)\n",
                 dyncolors_get_name_from_value(state_dynints->ints.a[j],
                                               dyncolors),
                 dyncolors_get_name_from_value(molecule->color_value,
                                               dyncolors),
                 dyncolors_get_name_from_value(component_color_value,
                                               dyncolors));
              return FALSE;
            }
            // State was not part of component initialization.
            // Warn about this, but it is not a fatal error.
            dynints_int_add(state_dynints->ints.a[j],
                            molecule->state_dynints_array.a[i]);
            DEBUG_PRINTF(0, "Adding state %s to %s(%s)\n",
                         dyncolors_get_name_from_value(state_dynints->ints.a[j],
                                                       dyncolors),
                         dyncolors_get_name_from_value(molecule->color_value,
                                                       dyncolors),
                         dyncolors_get_name_from_value(component_color_value,
                                                       dyncolors));
            error_printf("Adding state %s to %s(%s)\n",
                         dyncolors_get_name_from_value(state_dynints->ints.a[j],
                                                       dyncolors),
                         dyncolors_get_name_from_value(molecule->color_value,
                                                       dyncolors),
                         dyncolors_get_name_from_value(component_color_value,
                                                       dyncolors));
          }
        }
      }

      // Molecule definition is valid so far.
      return TRUE;
    }
  }

  DEBUG_PRINTF(0, "Unknown component %s(%s)\n",
               dyncolors_get_name_from_value(molecule->color_value,
                                             dyncolors),
               dyncolors_get_name_from_value(component_color_value,
                                             dyncolors));
  error_printf("Unknown component %s(%s)\n",
               dyncolors_get_name_from_value(molecule->color_value,
                                             dyncolors),
               dyncolors_get_name_from_value(component_color_value,
                                             dyncolors));
  // Component not found.
  return FALSE;
}

boolean
molecule_component_add(const int component_color_value,
                       const int component_multiplier,
                       dynints_t *dynint_states, dyncolors_t *dyncolors,
                       molecule_t *molecule)
{
  boolean status;
  size_t component_count;
  int *current_component_color_value;
  size_t i;

  // Search for existing component.
  component_count = molecule->component_count;
  for (i = 0,
       current_component_color_value = molecule->component_color_values.a;
       i < component_count;
       i++, current_component_color_value++)
  {
    if (component_color_value == *current_component_color_value)
    {
      if (0 == dynints_compare(dynint_states,
          molecule->state_dynints_array.a[i]))
      {
        DEBUG_PRINTF(0, "Duplicate components %s do not "
                     "have matching states\n",
                     dyncolors_get_name_from_value
                       (*current_component_color_value, dyncolors));
      }
      molecule->component_multipliers.a[i] += component_multiplier;
      return TRUE;
    }
  }

  // Component doesn't exist yet, so create new one.
  status = molecule_component_create(component_color_value,
                                     component_multiplier, dynint_states,
                                     molecule);
  if (FALSE == status)
  {
    return FALSE;
  }

  return TRUE;
}

boolean
molecule_component_create(const int component_color_value,
                          const int component_multiplier,
                          dynints_t *dynint_states,
                          molecule_t *molecule)
{
  boolean status;
  size_t index;

  // Need to create new component.
  status = molecule_component_alloc(molecule);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

  index = molecule->component_count - 1;

  molecule->component_color_values.a[index] = component_color_value;
  molecule->component_multipliers.a[index] = component_multiplier;
  if (NULL == molecule->state_dynints_array.a[index])
  {
    molecule->state_dynints_array.a[index] = dynints_create(0);
  }
  dynints_copy(dynint_states, molecule->state_dynints_array.a[index]);

  return TRUE;
}

boolean
molecule_component_alloc(molecule_t *molecule)
{
  int status;

  status = molecule_component_resize(molecule, molecule->component_count + 1);
  if (FALSE == status)
  {
    assert(0);
    return FALSE;
  }

  molecule->component_count++;

  return TRUE;
}

static boolean
molecule_component_resize(molecule_t *molecule, const size_t count)
{
  int status;

  status = dynarray_resize3(count, &molecule->component_alloc_count,
                            sizeof(int),
                            &molecule->component_color_values.v,
                            sizeof(int),
                            &molecule->component_multipliers.v,
                            sizeof(dynints_t *),
                            &molecule->state_dynints_array.v);
  if (0 == status)
  {
    return FALSE;
  }

  return TRUE;
}

void
molecule_destroy(molecule_t *molecule)
{
  int i;
  dynints_t **pcurrent_dynints;

  free(molecule->component_color_values.a);
  free(molecule->component_multipliers.a);
  for (i = 0, pcurrent_dynints = molecule->state_dynints_array.a;
       i < molecule->component_count;
       i++, pcurrent_dynints++)
  {
    if (NULL != *pcurrent_dynints)
    {
      free(*pcurrent_dynints);
    }
  }
  free(molecule->state_dynints_array.a);
  free(molecule);
}

