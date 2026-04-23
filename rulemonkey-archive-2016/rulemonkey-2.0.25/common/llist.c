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
** llist.c
**
** This library provides a simple linked list container.
*/

#include "llist.h"
#include <stdlib.h>
#include <assert.h>

/**
 * Return new linked list item initialized with arguments.
 */
llist_t *
llist_create(void *d, llist_t *next)
{
  llist_t *new;

  new = calloc(1, sizeof(llist_t));

  if (NULL == new)
  {
    // Memory allocation error.
    assert(0);
    return NULL;
  }

  new->d = d;
  new->next = next;

  return new;
}

/**
 * Create new linked list item at head of list, store p and update pllist.
 */
int
llist_push(void *p, llist_t **pllist)
{
  llist_t *new_llist;

  new_llist = llist_create(p, *pllist);
  if (NULL == new_llist)
  {
    assert(0);
    return 0;
  }

  *pllist = new_llist;
  return 1;
}

/**
 * Remove first item in linked list, returning contained pointer.
 */
void *
llist_shift(llist_t **pllist)
{
  llist_t *old_llist;
  void *data;

  old_llist = *pllist;

  *pllist = old_llist->next;

  data = old_llist->d;
  free(old_llist);

  return (data);
}

/**
 * Count number or items in linked list.
 */
size_t
llist_count(llist_t *llist)
{
  size_t count;

  count = 0;

  while(llist != NULL)
  {
    count++;
    llist = llist->next;
  }

  return count;
}

void
llist_destroy(llist_t *llist)
{
  if (NULL == llist)
  {
    // Nothing to do.
    return;
  }

  if (NULL != llist->next)
  {
    llist_destroy(llist->next);
  }

  free(llist);
}
