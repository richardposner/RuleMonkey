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
** molecule.h
**
** Defines a data structure for keeping track of a molecule.
*/
#ifndef MOLECULE_H
#define MOLECULE_H

typedef struct molecule_struct molecule_t;

#include "nauty.h"
#include "color.h"
#include "dynsymbols.h"
#include "dynints.h"
#include "dyncolors.h"
#include "dynarray.h"

struct molecule_struct {
  // Color of molecule.
  int color_value;

  // Dynamically allocated array of component colors.
  intarray_t component_color_values;

  // Dynamically allocated array of component counts.
  intarray_t component_multipliers;

  // Dynamically allocated array of state lists possible for each component;
  dynintsarray_t state_dynints_array;

  // Number of components currently in use.
  size_t component_count;

  // Number of components currently allocated.
  size_t component_alloc_count;
};

boolean molecule_setup(const int color_value,
                       const size_t initial_component_count,
                       molecule_t *new_molecule);
int *molecule_component_multipliers_table_create(molecule_t *molecule);
boolean molecule_component_verify
          (const int component_color_value,
           const int component_multiplier, 
           dynints_t *state_dynints,
           int *molecule_component_multipliers_table,
           dyncolors_t *dyncolors,
           molecule_t *molecule,
           int molecule_types_present);
boolean molecule_component_add(const int component_color_value,
                               const int component_multiplier,
                               dynints_t *dynint_states, 
                               dyncolors_t *dyncolors, molecule_t *molecule);
boolean molecule_component_create(const int component_color_value,
                                  const int component_multiplier, 
                                  dynints_t *state_dynints,
                                  molecule_t *molecule);
boolean molecule_component_alloc(molecule_t *molecule);
void molecule_destroy(molecule_t *molecule);


#endif
