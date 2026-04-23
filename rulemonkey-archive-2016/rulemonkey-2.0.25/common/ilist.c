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


#include "ilist.h"
#include "output.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

int
ilist_add(const int index, ilist_t **plist)
{
  ilist_t *new_ilist;

  new_ilist = ilist_create(index, *plist);
  if (NULL == new_ilist)
  {
    assert(0);
    return 0;
  }

  *plist = new_ilist;
  return 1;
}

ilist_t *
ilist_create(const int index, ilist_t *next)
{
  ilist_t *new_ilist;

  // Allocate new index.
  new_ilist = calloc(1, sizeof(ilist_t));
  if (NULL == new_ilist)
  {
    error_printf("Allocating new index");
    return NULL;
  }

  // Set structure members.
  new_ilist->index = index;
  new_ilist->next = next;

  return new_ilist;
}

int
ilist_pop(ilist_t **plist)
{
  ilist_t *old_ilist;
  int index;

  old_ilist = *plist;

  *plist = old_ilist->next;

  index = old_ilist->index;
  free(old_ilist);

  return(index);
}

/**
 * Recursively destroy all objects in list.
 */
void
ilist_destroy(ilist_t *ilist)
{
  if (NULL != ilist->next)
  {
    ilist_destroy(ilist->next);
  }

  free(ilist);
}
