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


#ifndef PARTICLESTOC_H
#define PARTICLESTOC_H

#include "world.h"
#include "dynps_species.h"
#include "output.h"

#define PARTICLESTOC_OUTPUT_ALL(ptype, level, world, log_time) \
  CHECKDEBUG(particlestoc_output_all, ptype, level, world, log_time)

boolean particlestoc_complete(world_t *world, dynps_species_t *dynps_species,
                              ps_particle_t **pps_particles,
                              size_t *pparticle_count);
void particlestoc_output_all(const output_type ptype, const int level,
                             world_t *world, const double log_time);
boolean particlestoc_run(world_t *world, const double logging_interval,
                         void (*report_callback)
                           (const output_type ptype, const int level,
                            world_t *, double),
                         const double max_time, const double volume,
                         const double time_delta, const double max_prob,
                         const double pseudo_factor,
                         const double pseudo_count, const double avogadro);

#endif
