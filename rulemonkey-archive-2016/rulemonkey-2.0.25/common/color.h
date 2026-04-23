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

#ifndef COLOR_H
#define COLOR_H

typedef struct color_struct color_t;
typedef enum color_type_enum color_type_t;

#define COLOR_WILDCARD -1

#define COLOR_ERROR_PREFIX "ERROR_"
#define COLOR_WILDCARD_PREFIX ""
#define COLOR_MOLECULE_PREFIX "molecule_"
#define COLOR_COMPONENT_PREFIX "component_"
#define COLOR_STATE_PREFIX "state_"
#define COLOR_COMPARTMENT_PREFIX "compartment_"

enum color_type_enum
{
  COLOR_TYPE_ERROR = 0,
  COLOR_TYPE_WILDCARD,
  COLOR_TYPE_MOLECULE,
  COLOR_TYPE_COMPONENT,
  COLOR_TYPE_STATE,
  COLOR_TYPE_COMPARTMENT
};

struct color_struct
{
  const char *name;

  // Name with appropriate prefix added.
  const char *normalized_name;

  int type;

  int value;
};

void color_setup(const char *name, color_type_t type,
                 const char *normalized_name,
                 const int value, color_t *new_color);
const char *color_normalize_name(const char *name, color_type_t type);
void color_destroy(color_t *color);

#endif
