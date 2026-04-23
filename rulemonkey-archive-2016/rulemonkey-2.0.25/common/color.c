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


#include <stdlib.h>
#include <string.h>
#include "color.h"

/*
** String pointed to by name is expected to have a lifetime
** equal to new_color, so it is not duplicated.
*/
void
color_setup(const char *name, color_type_t type,
            const char *normalized_name,
            const int value, color_t *new_color)
{
  new_color->name = name;
  new_color->type = type;
  new_color->normalized_name = normalized_name;
  new_color->value = value;
}

/*
** Returned string is allocated using malloc, so it should be cleaned up using
** free() when done with it.
*/
const char *
color_normalize_name(const char *name, color_type_t type)
{
  const char *prefix;
  size_t prefix_len;
  size_t name_len;
  size_t new_len;
  char *new_name;

  prefix_len = 0;
  prefix = "";
  name_len = strlen(name);

  switch (type)
  {
    case COLOR_TYPE_MOLECULE:
      prefix = COLOR_MOLECULE_PREFIX;
      // Don't count NUL character in prefix length.
      prefix_len = sizeof(COLOR_MOLECULE_PREFIX) - 1;
      break;
    case COLOR_TYPE_COMPONENT:
      prefix = COLOR_COMPONENT_PREFIX;
      // Don't count NUL character in prefix length.
      prefix_len = sizeof(COLOR_COMPONENT_PREFIX) - 1;
      break;
    case COLOR_TYPE_STATE:
      prefix = COLOR_STATE_PREFIX;
      // Don't count NUL character in prefix length.
      prefix_len = sizeof(COLOR_STATE_PREFIX) - 1;
      break;
    case COLOR_TYPE_COMPARTMENT:
      prefix = COLOR_COMPARTMENT_PREFIX;
      // Don't count NUL character in prefix length.
      prefix_len = sizeof(COLOR_COMPARTMENT_PREFIX) - 1;
      break;
    case COLOR_TYPE_WILDCARD:
      prefix = COLOR_WILDCARD_PREFIX;
      // Don't count NUL character in prefix length.
      prefix_len = sizeof(COLOR_WILDCARD_PREFIX) - 1;
      break;
    case COLOR_TYPE_ERROR:
      prefix = COLOR_ERROR_PREFIX;
      // Don't count NUL character in prefix length.
      prefix_len = sizeof(COLOR_ERROR_PREFIX) - 1;
      break;
  }

  // Add one for trailing character.
  new_len = name_len + prefix_len + 1;
  new_name = malloc(new_len);
  if (NULL == new_name)
  {
    return NULL;
  }

  // Begin new_name with prefix.
  memcpy(new_name, prefix, prefix_len);

  // Add name after prefix.
  memcpy(new_name + prefix_len, name, name_len);

  // Terminate string with NUL character.
  new_name[new_len - 1] = '\0';

  return new_name;
}


/* Frees any memory allocated inside color struct but not the struct itself */
void
color_destroy(color_t *color)
{
  free((char *)color->name);
}

