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
#include "dynmolecules.h"
#include "dynarray.h"
#include "constants.h"
#include "output.h"

void suppress_warnings24() {SUPPRESS_WARNINGS;}

// Local function declarations.
static int dynmolecules_resize(dynmolecules_t *dynmolecules,
                               const size_t count);

molecule_t *
dynmolecules_molecule_create(const int color_value,
                             const size_t component_count,
                             dynmolecules_t *dynmolecules)
{
  molecule_t *new_molecule;
  boolean status;

  // Make sure molecule does not exist yet.
  if (NULL != dynmolecules_molecule_find(color_value, dynmolecules))
  {
    // Molecule with given color already exists.
    error_printf("Molecule with color %d already exists.\n", color_value);
    return NULL;
  }

  new_molecule = dynmolecules_molecule_alloc(dynmolecules);
  if (NULL == new_molecule)
  {
    return NULL;
  }

  status = molecule_setup(color_value,
                          INITIAL_MOLECULE_COMPONENTS, new_molecule);
  if (FALSE == status)
  {
    dynmolecules->count--;
    return NULL;
  }

  return new_molecule;
}

molecule_t *
dynmolecules_molecule_alloc(dynmolecules_t *dynmolecules)
{
  int status;

  status = dynmolecules_resize(dynmolecules, dynmolecules->count + 1);
  if (0 == status)
  {
    assert(0);
    return NULL;
  }

  return &dynmolecules->molecules.a[dynmolecules->count++];
}

molecule_t *
dynmolecules_molecule_find(int color_value, dynmolecules_t *dynmolecules)
{
  int i;
  molecule_t *current_molecule;

  for (i = 0, current_molecule = dynmolecules->molecules.a;
       i < dynmolecules->count;
       i++, current_molecule++)
  {
    if (color_value == current_molecule->color_value)
    {
      // Found a match.
      return current_molecule;
    }
  }

  // Molecule was not found.
  return NULL;
}

dynmolecules_t *
dynmolecules_create(const size_t init_count)
{
  dynmolecules_t *new_dynmolecules;
  int status;

  new_dynmolecules = calloc(1, sizeof(dynmolecules_t));
  if (NULL == new_dynmolecules)
  {
    error_printf("Unable to allocate new dynmolecules.\n");
    return NULL;
  }

  status = dynarray_create(init_count, &new_dynmolecules->alloc_count,
                           sizeof(molecule_t),
                           &new_dynmolecules->molecules.v);
  if (0 == status)
  {
    // Error message output in dynarray_create.
    free(new_dynmolecules);
    return NULL;
  }

  new_dynmolecules->count = 0;

  return new_dynmolecules;
}

static int
dynmolecules_resize(dynmolecules_t *dynmolecules, const size_t count)
{
  int status;

  status = dynarray_resize(count, &dynmolecules->alloc_count,
                           sizeof(molecule_t),
                           &dynmolecules->molecules.v);
  if (0 == status)
  {
    // Error message output in dynarray_resize.
    return 0;
  }

  return 1;
}

void
dynmolecules_destroy(dynmolecules_t *dynmolecules)
{
  int i;
  molecule_t *current_molecule;

  for (i = 0, current_molecule = dynmolecules->molecules.a;
       i < dynmolecules->count;
       i++, current_molecule++)
  {
    molecule_destroy(current_molecule);
  }

  free(dynmolecules->molecules.a);
  free(dynmolecules);
}
