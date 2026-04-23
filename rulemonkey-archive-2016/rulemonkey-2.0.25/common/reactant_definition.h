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

#ifndef REACTANT_DEFINITION_H
#define REACTANT_DEFINITION_H

typedef struct reactant_definition_struct reactant_definition_t;

#include "dynrules.h"
#include "rule.h"

struct reactant_definition_struct
{
  // Linked list of rules that define reactant.
  rule_t *rule_list;

  // Name of rule to use.
  int rule_index;
};

reactant_definition_t *reactant_definition_create();

#endif
