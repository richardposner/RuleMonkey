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


#include "reactant_definition.h"
#include <assert.h>
#include "constants.h"

void suppress_warnings6() {SUPPRESS_WARNINGS;}

/*
** Returns pointer to newly allocated and initialized reactant definition.
*/
reactant_definition_t *
reactant_definition_create()
{
  reactant_definition_t *new_reactant_definition;

  new_reactant_definition = calloc(1, sizeof(reactant_definition_t));
  if (NULL == new_reactant_definition)
  {
    perror("Allocating new reactant definition");
    assert(0);
    return NULL;
  }

  return new_reactant_definition;
}

