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

#ifndef ILIST_H
#define ILIST_H

typedef struct ilist_struct ilist_t;
typedef union ilistarray_union ilistarray_t;

union ilistarray_union
{
  ilist_t *a;
  void *v;
};

struct ilist_struct
{
  int index;

  ilist_t *next;
};

int ilist_add(const int index, ilist_t **list);
ilist_t *ilist_create(const int index, ilist_t *next);
int ilist_pop(ilist_t **plist);
void ilist_destroy(ilist_t *ilist);

#endif
