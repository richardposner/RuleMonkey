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


#ifndef MYRAND_H
#define MYRAND_H

#include <dSFMT.h>
#include <stdlib.h>
#include <math.h>

#ifdef RANDOM_MT
#define MYRAND() myrand_mt()
//#define MYRAND_MAX(max) ((myrand_mt() * max) + 0.5)
#define MYRAND_INDEX(num) floor(myrand_mt() * num)
#define MYSRAND(seed) mysrand_mt(seed)
#else
//#define MYRAND() ran3()
#define MYRAND() (rand() / ((double)RAND_MAX + 1.0))
//#define MYRAND_INDEX(num) (ran3() * num)
#define MYRAND_INDEX(num) floor((MYRAND() * num))
//#define MYSRAND(seed) ran3_seed(seed)
#define MYSRAND(seed) srand(seed)
#endif

int mysrand_mt(const uint32_t seed);
double myrand_mt();

float ran3();
void ran3_seed(int seed);

#endif
