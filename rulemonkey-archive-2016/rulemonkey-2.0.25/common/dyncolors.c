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
#include "dyncolors.h"
#include "dynarray.h"
#include "output.h"

void suppress_warnings15() {SUPPRESS_WARNINGS;}

// Local function declarations.
static int dyncolors_check_size(dyncolors_t *dyncolors, const size_t count);

/*
** String pointed to by name must be allocated on heap, it will be
** freed using free() when needed.
*/
color_t *
dyncolors_color_create(const char *name, color_type_t type,
                       dyncolors_t *dyncolors)
{

  return real_dyncolors_color_create(name, type, NULL, dyncolors, NULL);
}

color_t *
real_dyncolors_color_create(const char *name, color_type_t orig_type,
                            dynspecies_t *dynspecies,
                            dyncolors_t *dyncolors, world_t *world)
{
  color_t *new_color;
  color_type_t type;

  if (0 == strcmp("*", name))
  {
    type = COLOR_TYPE_WILDCARD;
  }
  else {
    type = orig_type;
  }
  new_color = dyncolors_color_create_any(name, type, dynspecies, dyncolors,
                                         world);

  return new_color;
}

/*
** String pointed to by name must be allocated on heap, it will be
** freed using free() when needed.
*/
color_t *
dyncolors_color_create_any(const char *orig_name, color_type_t type,
                           dynspecies_t *dynspecies, dyncolors_t *dyncolors,
                           world_t *world)
{
  color_t *new_color;
  const char *normalized_name;

  new_color = dyncolors_color_name_type_lookup(orig_name, type, dyncolors);
  if (NULL != new_color)
  {
    // Color already exists.
    return new_color;
  }

  // Need to create new color.
  new_color = dyncolors_color_alloc(dyncolors);
  if (NULL == new_color)
  {
    return NULL;
  }

  // Return value normalized_name is allocated on the heap, so
  // need to be careful to make sure it will get cleaned up.
  normalized_name = color_normalize_name(orig_name, type);
  if (NULL == normalized_name)
  {
    return NULL;
  }

  // Color_setup uses normalized_name as well as name, so don't free it.
  color_setup(strdup(orig_name), type, normalized_name,
              (int)dyncolors->count - 1, new_color);

  return new_color;
}

void
dyncolors_output(const output_type ptype, const int level,
                 dyncolors_t *dyncolors, const char *prefix,
                 const char *postfix, boolean value, boolean prompt)
{
  int i;
  color_t *current_color;

  CUSTOM_PRINTF(ptype, level, "Colors:\n");

  for (i = 0, current_color = dyncolors->colors.a;
       i < dyncolors->count;
       i++, current_color++)
  {
    if (TRUE == value)
      {
        // Need to include integer value of color.
        CUSTOM_PRINTF(ptype, level, "  %s%s%s = %d\n", prefix,
                      current_color->normalized_name, postfix,
                      current_color->value);
      }
    else
      {
        CUSTOM_PRINTF(ptype, level, "  %s%s%s\n",
                      prefix, current_color->normalized_name, postfix);
      }
  }
}

color_t *
dyncolors_color_alloc(dyncolors_t *dyncolors)
{
  int status;

  status = dyncolors_check_size(dyncolors, dyncolors->count + 1);
  if (FALSE == status)
  {
    assert(0);
    return NULL;
  }

  return &dyncolors->colors.a[dyncolors->count++];
}

color_t *
dyncolors_color_normalized_name_lookup(const char *normalized_name,
                                       dyncolors_t *dyncolors)
{
  int i;
  color_t *current_color;

  if (NULL == normalized_name)
  {
    assert(0 && "dyncolors_color_lookup got an empty normalized_name pointer");
    return NULL;
  }

  current_color = dyncolors->colors.a;

  for (i = 0, current_color = dyncolors->colors.a;
       i <dyncolors->count;
       i++, current_color++)
  {
    if (0 == strcmp(normalized_name, current_color->normalized_name))
    {
      // Found a match.
      return current_color;
    }
  }

  // Name was not found.
  return NULL;
}

color_t *
dyncolors_color_name_type_lookup(const char *name, color_type_t type,
                                 dyncolors_t *dyncolors)
{
  int i;
  color_t *current_color;

  if (NULL == name)
  {
    assert(0 && "dyncolors_color_lookup got an empty name pointer");
    return NULL;
  }

  for (i = 0, current_color = dyncolors->colors.a;
       i <dyncolors->count;
       i++, current_color++)
  {
    if ((type == current_color->type) &&
        (0 == strcmp(name, current_color->name)))
    {
      // Found a match.
      return current_color;
    }
  }

  // Name was not found.
  return NULL;
}

color_t *
dyncolors_color_value_lookup(const int value, dyncolors_t *dyncolors)
{
  int i;
  color_t *current_color;

  for (i = 0, current_color = dyncolors->colors.a;
       i <dyncolors->count;
       i++, current_color++)
  {
    if (value == current_color->value)
    {
      // Found a match.
      return current_color;
    }
  }

  // Value was not found.
  return NULL;
}

color_type_t
dyncolors_get_type_from_value(const int value, dyncolors_t *dyncolors)
{
  color_t *current_color;

  current_color = dyncolors_color_value_lookup(value, dyncolors);
  if (NULL == current_color)
  {
    return COLOR_TYPE_ERROR;
  }

  return current_color->type;
}

const char *
dyncolors_get_name_from_value(const int value, dyncolors_t *dyncolors)
{
  color_t *current_color;

  current_color = dyncolors_color_value_lookup(value, dyncolors);
  if (NULL == current_color)
  {
    return NULL;
  }

  return current_color->name;
}

dyncolors_t *
dyncolors_create(const size_t init_count)
{
  dyncolors_t *new_dyncolors;
  int status;
  color_t *color;

  // All fields need to be initialized to zero anyways, so use calloc.
  new_dyncolors = calloc(1, sizeof(dyncolors_t));
  if (NULL == new_dyncolors)
  {
    error_printf("Unable to allocate new dyncolors.\n");
    return NULL;
  }

  status = dynarray_create(init_count, &new_dyncolors->alloc_count,
                           sizeof(color_t), &new_dyncolors->colors.v);
  if (0 == status)
  {
    // Error message output in dynarray_create.
    error_printf("Unable to allocate dyncolors->colors.\n");
    free(new_dyncolors);
    return NULL;
  }

  // COLOR_WILDCARD always exists.
  color = dyncolors_color_create_any("*", COLOR_TYPE_WILDCARD,
                                     NULL, new_dyncolors, NULL);
  if (NULL == color)
  {
    // Unable to create wildcard color.
    error_printf("Unable to create wildcard color.\n");
    free(new_dyncolors->colors.a);
    free(new_dyncolors);
    return NULL;
  }
  new_dyncolors->colors.a[0].value = COLOR_WILDCARD;

  return new_dyncolors;
}

static int
dyncolors_check_size(dyncolors_t *dyncolors, const size_t count)
{
  int status;

  status = dynarray_resize(count, &dyncolors->alloc_count, sizeof(color_t),
                           &dyncolors->colors.v);
  if (0 == status)
  {
    // Error message output in dynarray_resize.
    return 0;
  }

  return 1;
}

void
dyncolors_destroy(dyncolors_t *dyncolors)
{
  int i;
  color_t *current_color;

  for (i = 0, current_color = dyncolors->colors.a;
       i <dyncolors->count;
       i++, current_color++)
  {
    color_destroy(current_color);
  }

  free(dyncolors->colors.a);
  free(dyncolors);
}
