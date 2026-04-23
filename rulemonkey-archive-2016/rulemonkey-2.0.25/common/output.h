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

#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdio.h>
#include <math.h>

enum output_type_t
{
  OUTPUT_STANDARD,
  OUTPUT_BNG_STANDARD,
  OUTPUT_DEBUG,
  OUTPUT_ERROR
};

typedef enum output_type_t output_type;

extern int current_debug_level;
extern FILE *my_out;
extern FILE *my_bng_out;
extern FILE *my_dbg;
extern FILE *my_err;

#define CUSTOM_OUTPUT_CHECK(level) ((current_debug_level >= level) ? 1 : 0)

#define CHECKDEBUG(func, ptype, level, ...) \
   if (current_debug_level >= level) func(ptype, level, __VA_ARGS__)

#define DEBUG_PRINTF(level, ...) \
   CHECKDEBUG(custom_printf, OUTPUT_DEBUG, level, __VA_ARGS__)

#define CUSTOM_PRINTF(ptype, level, ...) \
   CHECKDEBUG(custom_printf, ptype, level, __VA_ARGS__)

int output_begin(const char *raw_name);
void output_end();
void set_debug_level(const int debug_level);

FILE *custom_get_fh(const output_type ptype);
void custom_printf(const output_type ptype, const int level,
                   const char *fmt, ...)
       __attribute__ ((format (printf, 3, 4)));
void bng_standard_printf(const char *fmt, ...)
       __attribute__ ((format (printf, 1, 2)));
void standard_printf(const char *fmt, ...)
       __attribute__ ((format (printf, 1, 2)));
void debug_printf(const int debug_level, const char *fmt, ...)
       __attribute__ ((format (printf, 2, 3)));
void error_printf(const char *fmt, ...)
       __attribute__ ((format (printf, 1, 2)));

void output_procinfo(const output_type ptype, const int level);

#endif /* OUTPUT_H */
