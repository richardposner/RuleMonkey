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


#ifndef PS_REACTION_DEFINITION_H
#define PS_REACTION_DEFINITION_H

#include "reaction_definition.h"

void ps_reaction_definition_init_all_prob
  (const double particle_count, const double pseudo_particle_count,
   double *time_delta, const double volume,
   reaction_definition_t *reaction_definition_list, const double max_prob,
   const double avogadro);

#endif
