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

/*
** dynmolecules.h
**
** Data structure for keeping tack of molecules.
*/
#ifndef DYNMOLECULES_H
#define DYNMOLECULES_H

typedef struct dynmolecules_struct dynmolecules_t;

typedef union moleculearray_union moleculearray_t;

#include "molecule.h"

union moleculearray_union
{
  molecule_t *a;
  void *v;
};

struct dynmolecules_struct {
  // Dynamically allocated array of molecules.
  moleculearray_t molecules;

  // Number of molecules currently in use.
  size_t count;

  // Number of molecules currently allocated.
  size_t alloc_count;
};

molecule_t *dynmolecules_molecule_create(const int color_value,
                                         const size_t component_count,
                                         dynmolecules_t *dynmolecules);
molecule_t *dynmolecules_molecule_alloc(dynmolecules_t *dynmolecules);
molecule_t *dynmolecules_molecule_find(int color_value,
                                       dynmolecules_t *dynmolecules);
dynmolecules_t *dynmolecules_create(const size_t init_count);
void dynmolecules_destroy(dynmolecules_t *dynmolecules);

#endif
