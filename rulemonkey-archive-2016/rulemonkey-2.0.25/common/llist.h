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

/*
** llist.h
**
** This library provides a simple linked list container.
*/
#ifndef LLIST_H
#define LLIST_H

typedef struct llist_struct llist_t;

#include "llist.h"
#include <stdio.h>

struct llist_struct
{
  // Pointer to object being contained.
  void *d;

  // Next in list.
  llist_t *next;
};

llist_t *llist_create(void *d, llist_t *next);
int llist_push(void *p, llist_t **pllist);
void *llist_shift(llist_t **pllist);
void llist_destroy(llist_t *llist);
size_t llist_count(llist_t *llist);

#endif
