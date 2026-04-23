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

#include "parser.h"
#include "parser_lib.h"
#include "dyncolors.h"
#include "output.h"
#include "string.h"


void suppress_warnings32() {SUPPRESS_WARNINGS;}

extern FILE *yyin;

int main(int argc, char **argv)
{
  int status;
  char *filename;

  //yydebug = 1;

  if (argc > 1)
  {
    FILE *file;

    filename = argv[1];

    if (0 == strncasecmp(filename, "-v", 2))
    {
      printf(PACKAGE_STRING "\n");
      exit(0);
    }

    file = fopen(filename, "r");
    if (!file)
    {
      fprintf(stderr, "Could not open %s\n", argv[1]);
      exit(1);
    }

    yyin = file;
  }
  else
  {
    filename = NULL;
  }

  status = output_begin(filename);
  if (0 == status)
  {
    return (1);
  }

  status = setup_parse_env();
  if (FALSE == status)
  {
    return (1);
  }

  yyparse();

  output_procinfo(OUTPUT_DEBUG, 0);

  output_end();

  return (0);
}

