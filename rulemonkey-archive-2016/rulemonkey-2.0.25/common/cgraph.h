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

/**
** @file cgraph.h
**
** @brief This library provides a convenient way for using colored graphs with
** the NAUTY library.
**
** @author Joshua Colvin (Translation Genomics Research Institute)
**
** Copyright (C) 2007 Joshua Colvin, Northern Arizona University,
** Translational Genomics Research Institute. All rights reserved.
*/

#ifndef CGRAPH_H
#define CGRAPH_H

typedef struct cgraph_struct cgraph_t;
typedef union grapharray_union grapharray_t;

#include "../nauty22/nauty.h"
#include "dynints.h"
#include "dynarray.h"
#include "dyncolors.h"
#include "output.h"

#define CGRAPH_OUTPUT_BNGL(ptype, level, cg, dyncolors) \
   CHECKDEBUG(cgraph_output_bngl, ptype, level, cg, dyncolors)

#define CGRAPH_OUTPUT(ptype, level, cg, dyncolors) \
   CHECKDEBUG(cgraph_output, ptype, level, cg, dyncolors)

// To prevent type-punning warnings.
union grapharray_union
{
  graph *a;
  void *v;
};

struct cgraph_struct
{
  int n; // Number of vertices.
  int n_alloc; // Allocated number of vertices, n_alloc >= n.
  int m; // Number of setwords in sets.

  // Canonically labelled graph.
  grapharray_t g;

  // Canonical labelling of the graph
  intarray_t lab;
  intarray_t ptn;

  // Color value: lab[i] => colors[i]
  intarray_t colors;

  // State(s) assiciated with each node: lab[i] => state_lists[i]
  // Individual states are identified by a 'color'.
  dynintsarray_t state_dynints_array;

  // Used by species to enforce exclusive adjacency.
  intarray_t adjacent_counts;

  // Orbits
  intarray_t orbits;
};

cgraph_t *cgraph_create(int n);
void cgraph_output_bngl(const output_type ptype, const int level,
                   cgraph_t *cg, dyncolors_t *dyncolors);
void cgraph_output(const output_type ptype, const int level,
                   cgraph_t *cg, dyncolors_t *dyncolors);
cgraph_t *cgraph_create_raw(int n);
boolean cgraph_realloc(cgraph_t *a, int n);
void cgraph_reset(cgraph_t *a);
void cgraph_fix_ptn(cgraph_t *a);
int cgraph_normalize_partitions(cgraph_t *a);
boolean cgraph_normalize(cgraph_t *a);
boolean cgraph_free(cgraph_t *a);
boolean cgraph_destroy(cgraph_t *a);
cgraph_t *cgraph_copy(cgraph_t *src);
boolean cgraph_ncopy(cgraph_t *src, cgraph_t *dest);
boolean cgraph_canon (cgraph_t *a);
boolean cgraph_compare_nodes (cgraph_t *a, int e, int f);
boolean cgraph_compare (cgraph_t *a, cgraph_t *b);
boolean cgraph_link (cgraph_t *a, int v_a, cgraph_t *b, int v_b,
                 cgraph_t *c);
boolean cgraph_append (cgraph_t *a, int *a_renum,
                       cgraph_t *b, int *b_renum,
                       cgraph_t *c);
boolean cgraph_unlink (cgraph_t *a, int v_a, int v_b,
                   cgraph_t *c, cgraph_t *d);
boolean cgraph_seperate (cgraph_t *a, cgraph_t *c, cgraph_t *d);
boolean cgraph_seperate_raw (cgraph_t *a, cgraph_t *c, cgraph_t *d,
                             int *c_renum, int *d_renum);
int cgraph_add_node(cgraph_t *a, int color);
void cgraph_connect_nodes(cgraph_t *a, int e, int f);



#endif
