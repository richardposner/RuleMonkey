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


#ifndef REPORT_DEFINITION_CONTAINER_H
#define REPORT_DEFINITION_CONTAINER_H

#include "report_definition.h"

typedef struct report_definition_container_struct report_definition_container_t;

struct report_definition_container_struct
{
  // Report definition that is associated with this container.
  report_definition_t *report_definition;

  // Multiplier involved with this container.
  int multiplier;
};

#endif
