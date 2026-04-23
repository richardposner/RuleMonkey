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
** @file cgraph.c
**
** @brief This library provides a convenient way for using colored graphs with
** the NAUTY library.
**
** @author Joshua Colvin (Translation Genomics Research Institute)
**
** Copyright (C) 2008 Northern Arizona University and
** The Translational Genomics Research Institute. All rights reserved.
*/

#include "cgraph.h"
#include "../nauty22/naututil.h"
#include "assert.h"
#include "constants.h"
#include "dynarray.h"
#include "output.h"

void suppress_warnings3() {SUPPRESS_WARNINGS;}

static int cgraph_get_m(int n);
static void cgraph_reset_raw(cgraph_t *a);
static boolean
cgraph_seperate_raw_raw (cgraph_t *a, cgraph_t *c, cgraph_t *d,
                         set *work_set, set *c_set, int *c_renum, int *d_renum);
static boolean cgraph_realloc_raw(cgraph_t *a, int n);

/**
** This function creates an empty graph structure.
**
** @param n initial size of graph
** @return cgraph object.
*/
cgraph_t *
cgraph_create(int n)
{
  cgraph_t *a;

  a = cgraph_create_raw(n);

  // Initialize memory to all zeros.
  if (NULL != a)
  {
    cgraph_reset(a);
  }

  return a;
}

/**
** Output the current state of the cgraph object in human-readable form.
**
** @param ptype Where to send output to
** @param level Debug level
** @param cg Object to output
** @param dyncolors Used to map label index to name.
*/
void
cgraph_output_bngl(const output_type ptype, const int level,
                   cgraph_t *cg, dyncolors_t *dyncolors)
{
  dynints_t *connections;
  graph *g;
  int n;
  int m;
  int i;
  int molecule_count;

  n = cg->n;
  m = cg->m;
  g = cg->g.a;

  connections = dynints_create(10);
  molecule_count = 0;

  // Go through each node.
  for (i = 0; i < n; i++)
  {
    color_t *molecule_color;
    int component_count;
    int j;
    set *molecule_gv;

    molecule_color = dyncolors_color_value_lookup(cg->colors.a[i], dyncolors);

    if (COLOR_TYPE_MOLECULE != molecule_color->type)
    {
      // Nothing to do here.
      continue;
    }

    // Found a molecule, so start by printing the name.
    if (0 != molecule_count)
    {
      CUSTOM_PRINTF(ptype, level, ".");
    }
    molecule_count++;
    CUSTOM_PRINTF(ptype, level, "%s(", molecule_color->name);

    molecule_gv = GRAPHROW(g, i, m);

    component_count = 0;

    // Go through each connection molecule has.
    for (j = -1; (j = nextelement(molecule_gv, m, j)) >= 0;)
    {
      color_t *component_color;
      dynints_t *state_dynints_array;
      int *ints;
      size_t count;
      int k;
      set *component_gv;

      component_color =
        dyncolors_color_value_lookup(cg->colors.a[j], dyncolors);

      assert(COLOR_TYPE_COMPONENT == component_color->type);

      // If not first component, need to output comma.
      if (0 != component_count)
      {
        CUSTOM_PRINTF(ptype, level, ",");
      }
      component_count++;

      // Output component name.
      CUSTOM_PRINTF(ptype, level, "%s", component_color->name);

      component_gv = GRAPHROW(g, j, m);

      // Output states.
      state_dynints_array = cg->state_dynints_array.a[j];
      ints = state_dynints_array->ints.a;
      count = state_dynints_array->count;
      for (k = 0; k < count; k++)
      {
        color_t *state_color;

        state_color = dyncolors_color_value_lookup(ints[k], dyncolors);

        if (NULL == state_color)
        {
          CUSTOM_PRINTF(ptype, level, "~%d*ERROR*", ints[k]);
        }
        else
        {
          CUSTOM_PRINTF(ptype, level, "~%s", state_color->name);
        }
      }

      // Go through each connection component has.
      for (k = -1; (k = nextelement(component_gv, m, k)) >= 0;)
      {
        color_t *connection_color;
        int connection_count;
        int smallest;

        // Go through each connection component has.
        if (k == i)
        {
          // Skip connection back to molecule.
          continue;
        }

        connection_color =
          dyncolors_color_value_lookup(cg->colors.a[k], dyncolors);

        assert(COLOR_TYPE_COMPONENT == connection_color->type);

        if (k > j)
        {
          smallest = j;
        }
        else
        {
          smallest = k;
        }

        connection_count = dynints_int_find(smallest, connections);

        if (-1 == connection_count)
        {
          // Connection not encountered yet.
          // Save current node so that we can lookup the connection later.
          dynints_int_add(smallest, connections);

          connection_count = dynints_int_find(smallest, connections);
        }

        // Output connection.
        // connection count is 0 based, make connections 1 based.
        CUSTOM_PRINTF(ptype, level, "!%d", connection_count + 1);
      }
    }

    // Finished current molecule.
    CUSTOM_PRINTF(ptype, level, ")");
  }


  dynints_destroy(connections);
}

/**
** Output the current state of the cgraph object in human-readable form.
**
** @param ptype Where to send output to
** @param level Debug level
** @param cg Object to output
** @param dyncolors Used to map label index to name.
*/
void
cgraph_output(const output_type ptype, const int level,
              cgraph_t *cg, dyncolors_t *dyncolors)
{
  int orig_labelorg;

  orig_labelorg = labelorg;
#if GRAPH_BASE_ONE
  labelorg = 1;
#endif

  int i;
  color_t *color;

  putgraph(custom_get_fh(ptype), cg->g.a, 80, cg->m, cg->n);

  CUSTOM_PRINTF(ptype, level, "Colors:\n");
  for (i = 0; i < cg->n; i++)
  {
    color = dyncolors_color_value_lookup(cg->colors.a[i], dyncolors);
    if (NULL != color)
    {
      CUSTOM_PRINTF(ptype, level, "  %d : %s\n",
                    i + labelorg, color->normalized_name);
    }
    else
    {
      // Unable to find color name.
      CUSTOM_PRINTF(ptype, level, "  %d : %d\n", i + labelorg, cg->colors.a[i]);
    }
  }

  CUSTOM_PRINTF(ptype, level, "States:\n");
  for (i = 0; i < cg->n; i++)
  {
    dynints_t *state_dynints_array;
    int *ints;
    size_t count;
    size_t j;

    CUSTOM_PRINTF(ptype, level, "  %d : ", i + labelorg);

    state_dynints_array = cg->state_dynints_array.a[i];
    ints = state_dynints_array->ints.a;
    count = state_dynints_array->count;

    for (j = 0; j < count; j++)
    {
      color = dyncolors_color_value_lookup(ints[j], dyncolors);
      if (NULL != color)
      {
        CUSTOM_PRINTF(ptype, level, " %s", color->normalized_name);
      }
      else
      {
        // Unable to find color name.
        CUSTOM_PRINTF(ptype, level, "%d ", ints[j]);
      }
    }
    CUSTOM_PRINTF(ptype, level, "\n");
  }

  labelorg = orig_labelorg;
}

/**
** Calculate size of m needed to store graph with the given number of vertices.
**
** @param n number of vertices in graph
** @return computed value for m
*/
static int
cgraph_get_m(int n)
{
  return (n + WORDSIZE - 1) / WORDSIZE;
}

/**
** Create a new colored graph, do not reset memory to zero.
**
** @param n number of vertices to initialize graph with
** @return newly created cgraph object, not initialized
*/
cgraph_t *
cgraph_create_raw(int n)
{
  cgraph_t *a;
  boolean status;
  int i;
  dynints_t **current_state_dynints;

  a = calloc(1, sizeof(cgraph_t));
  if (NULL == a)
  {
    assert(0);
    return NULL;
  }

  if (n < 1)
  {
    // Make sure n is at least 1 so that memory block is allocated.
    n = 1;
  }

  a->n = n;
  a->n_alloc = n;
  a->m = cgraph_get_m(a->n);

  status = dynarray_realloc(0, (size_t)(a->m * a->n),
                            sizeof(graph), &a->g.v);
  if (FALSE == status)
  {
    assert(0);
    free(a);
    return NULL;
  }
  status = dynarray_realloc(0, (size_t)(a->n), sizeof(int), &a->lab.v);
  if (FALSE == status)
  {
    assert(0);
    free(a->g.v);
    a->g.v = NULL;
    free(a);
    return NULL;
  }
  status = dynarray_realloc(0, (size_t)(a->n), sizeof(int), &a->ptn.v);
  if (FALSE == status)
  {
    assert(0);
    free(a->lab.a);
    a->lab.a = NULL;
    free(a->g.a);
    a->g.a = NULL;
    free(a);
    return NULL;
  }
  status = dynarray_realloc(0, (size_t)(a->n),
                            sizeof(int), &a->colors.v);
  if (FALSE == status)
  {
    assert(0);
    free(a->ptn.a);
    a->ptn.a = NULL;
    free(a->lab.a);
    a->lab.a = NULL;
    free(a->g.a);
    a->g.a = NULL;
    free(a);
    return NULL;
  }
  status = dynarray_realloc(0, (size_t)(a->n),
                            sizeof(int), &a->adjacent_counts.v);
  if (FALSE == status)
  {
    assert(0);
    free(a->colors.a);
    a->colors.a = NULL;
    free(a->ptn.a);
    a->ptn.a = NULL;
    free(a->lab.a);
    a->lab.a = NULL;
    free(a->g.a);
    a->g.a = NULL;
    free(a);
    return NULL;
  }
  status = dynarray_realloc(0, (size_t)(a->n),
                            sizeof(int), &a->orbits.v);
  if (FALSE == status)
  {
    assert(0);
    free(a->adjacent_counts.a);
    a->adjacent_counts.a = NULL;
    free(a->colors.a);
    a->colors.a = NULL;
    free(a->ptn.a);
    a->ptn.a = NULL;
    free(a->lab.a);
    a->lab.a = NULL;
    free(a->g.a);
    a->g.a = NULL;
    free(a);
    return NULL;
  }
  status = dynarray_realloc(0, (size_t)(a->n),
                            sizeof(dynints_t *),
                            &a->state_dynints_array.v);
  if (FALSE == status)
  {
    assert(0);
    free(a->orbits.a);
    free(a->adjacent_counts.a);
    a->adjacent_counts.a = NULL;
    free(a->colors.a);
    a->colors.a = NULL;
    free(a->ptn.a);
    a->ptn.a = NULL;
    free(a->lab.a);
    a->lab.a = NULL;
    free(a->g.a);
    a->g.a = NULL;
    free(a);
    return NULL;
  }
  for (i = 0, current_state_dynints = a->state_dynints_array.a;
       i < a->n_alloc;
       i++, current_state_dynints++)
  {
    *current_state_dynints = dynints_create(INITIAL_DYNINTS);
    if (NULL == *current_state_dynints)
    {
      int j;

      // Memory allocation failed, free everything we've done so far.
      for (j = i - 1, current_state_dynints--;
           j >= 0;
           j--)
      {
        dynints_destroy(*current_state_dynints);
      }
      status = FALSE;
      break;
    }
  }
  if (FALSE == status)
  {
    assert(0);
    free(a->state_dynints_array.a);
    free(a->orbits.a);
    free(a->adjacent_counts.a);
    a->adjacent_counts.a = NULL;
    free(a->colors.a);
    a->colors.a = NULL;
    free(a->ptn.a);
    a->ptn.a = NULL;
    free(a->lab.a);
    a->lab.a = NULL;
    free(a->g.a);
    a->g.a = NULL;
    free(a);
    return NULL;
  }

  return a;
}

/**
** Reset cgraph and make sure cgraph is large enough to hold graph.
**
** @param a cgraph object to reset
** @param n size of new graph
** @return 0 if error, 1 otherwise.
*/
boolean
cgraph_realloc_clear(cgraph_t *a, int n)
{
  boolean status;

  cgraph_reset(a);

  status = cgraph_realloc_raw(a, n);
  if (FALSE == status)
  {
    return FALSE;
  }

  cgraph_reset_raw(a);

  return TRUE;
}

/**
** Make sure cgraph is large enough to hold graph.
**
** @param a cgraph object to check
** @param n number of vertices graph needs to hold
** @return 0 if error, 1 otherwise
*/
boolean
cgraph_realloc(cgraph_t *a, int n)
{
  boolean status;
  int old_m;
  int old_n;

  old_m = a->m;
  old_n = a->n;

  status = cgraph_realloc_raw(a, n);
  if (FALSE == status) return FALSE;

  if (old_m == a->m)
  {
    // Nothing left to do.
    return TRUE;
  }

  // The value of m changed, which means each node takes up a different
  // amount of memory.
  if (old_m > a->m)
  {
    // Haven't implemented shrinking graph yet.
    assert(0);
    return FALSE;
  }
  else /* old_m < a->m */
  {
    graph *old_node;
    graph *new_node;
    int i;

    // Nodes grew larger.  Reposition each node, starting at the last node
    // so that we don't have to have two copies of graph.
    // Use old_n because there is no need to copy new empty nodes.
    for (i = 0, old_node = a->g.a + (old_n * old_m),
                new_node = a->g.a + (old_n * a->m);
         i < a->n;
         i++, old_node -= old_m,
              new_node -= a->m)
    {
      // Copy old node to new location.
      // Use memmove because regions may overlap.
      memmove(new_node, old_node, sizeof(graph) * old_m);

      // Clear new node memory that wasn't initialized by old node memory.
      memset(new_node + old_m, 0, sizeof(graph) * (a->m - old_m));
    }
  }

  return TRUE;
}

/**
** Resize graph if needed, do not reset memory to zero.
**
** @param a cgraph object to check
** @param n number of vertices graph needs to hold
** @return 0 if error, 1 otherwise
*/
static boolean
cgraph_realloc_raw(cgraph_t *a, int n)
{
  size_t n_alloc;
  size_t m_alloc;
  size_t old_m_alloc;
  size_t old_n_alloc;
  boolean status;
  int i;
  dynints_t **current_state_dynints;

  if (0 == n)
  {
    a->n = 0;
    a->m = 0;
    return TRUE;
  }

  a->n = n;
  a->m = cgraph_get_m(a->n);

  old_n_alloc = a->n_alloc;
  old_m_alloc = (old_n_alloc + WORDSIZE - 1) / WORDSIZE;

  n_alloc = dynarray_check_size((size_t)n, (size_t)a->n_alloc);

  if (n_alloc > 0)
  {
    // Need to allocate more memory.
    a->n_alloc = n_alloc;
    m_alloc = (n_alloc + WORDSIZE - 1) / WORDSIZE;

    status = dynarray_realloc((size_t)(old_m_alloc * old_n_alloc),
                              (size_t)(m_alloc * n_alloc),
                              sizeof(graph), &a->g.v);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }
    status = dynarray_realloc((size_t)(old_n_alloc), (size_t)(n_alloc),
                              sizeof(int), &a->lab.v);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }
    status = dynarray_realloc((size_t)(old_n_alloc), (size_t)(n_alloc),
                              sizeof(int), &a->ptn.v);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }
    status = dynarray_realloc((size_t)(old_n_alloc), (size_t)(n_alloc),
                              sizeof(int), &a->colors.v);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }
    status = dynarray_realloc((size_t)(old_n_alloc), (size_t)(n_alloc),
                              sizeof(int), &a->adjacent_counts.v);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }
    status = dynarray_realloc((size_t)(old_n_alloc), (size_t)(n_alloc),
                              sizeof(int), &a->orbits.v);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }
    status = dynarray_realloc((size_t)(old_n_alloc), (size_t)(n_alloc),
                              sizeof(dynints_t *),
                              &a->state_dynints_array.v);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }
    for (i = old_n_alloc,
         current_state_dynints = a->state_dynints_array.a + old_n_alloc;
         i < n_alloc;
         i++, current_state_dynints++)
    {
      *current_state_dynints = dynints_create(INITIAL_DYNINTS);
      if (NULL == *current_state_dynints)
      {
        return FALSE;
      }
    }

  } 
  return TRUE;
}

/**
** Zero out all memory related to cgraph.
**
** @param a object to manipulate
*/
void
cgraph_reset(cgraph_t *a)
{
  cgraph_reset_raw(a);

  a->m = 0;
  a->n = 0;
}

/**
** Erase all graph info, but keep graph size intact.
**
** @param a object to manipulate
*/
static void
cgraph_reset_raw(cgraph_t *a)
{
  int m_alloc;
  int i;
  dynints_t **pcurrent_states;

  m_alloc = cgraph_get_m(a->n_alloc);

  memset(a->g.a,      0, (size_t)(m_alloc * a->n_alloc) * sizeof(graph));
  memset(a->lab.a,    0, (size_t)a->n_alloc * sizeof(int));
  memset(a->ptn.a,    0, (size_t)a->n_alloc * sizeof(int));
  memset(a->colors.a, 0, (size_t)a->n_alloc * sizeof(int));
  memset(a->adjacent_counts.a, 0, (size_t)a->n_alloc * sizeof(int));
  memset(a->orbits.a, 0, (size_t)a->n_alloc * sizeof(int));
  for (i = 0, pcurrent_states = a->state_dynints_array.a;
       i < a->n_alloc;
       i++, pcurrent_states++)
  {
    (*pcurrent_states)->count = 0;
  }
}

/**
** Add node with given color to cgraph.
**
** Note: Need to normalize the graph partitions once graph is
**       fully built.
**
** @param a object to manipulate
** @param color 
** @return index of new node.
*/
int
cgraph_add_node(cgraph_t *a, int color)
{
  boolean status;
  int new_index;

  // Make sure enough memory allocated for additional node.
  status = cgraph_realloc(a, a->n + 1);
  if (FALSE == status)
  {
    assert(0);
    return 0;
  }

  new_index = a->n - 1;

  // Set color value for new node.
  a->colors.a[new_index] = color;

  return new_index;
}

/**
** Connect two nodes together.
**
** @param a object to manipulate
** @param e first node to connect
** @param f second node to connect
*/
void
cgraph_connect_nodes(cgraph_t *a, int e, int f)
{
  set *e_gv;
  set *f_gv;

#ifdef AGGRESSIVE_CHECK
assert(e <= a->n);
assert(f <= a->n);
#endif
  e_gv = GRAPHROW(a->g.a, e, a->m);
  f_gv = GRAPHROW(a->g.a, f, a->m);

  ADDELEMENT(e_gv, f);
  ADDELEMENT(f_gv, e);
}

/**
** Create partitions such that colors are normalized (in numerical order).
** Example: 2010602626 becomes 0001222666
**
** @param a object to manipulate
*/
int
cgraph_normalize_partitions(cgraph_t *a)
{
  int node_count;
  int current_color;
  int label_index;
  int node_index;

  // Group all colors in contiguous blocks
  node_count = a->n;
  current_color = 0;
  label_index = 0;

  // Check that there are no 0s.
  for (node_index = 0;
       node_index < node_count;
       node_index++)
  {
    if (0 == a->colors.a[node_index])
    {
      error_printf("Invalid color found in cgraph.\n");
      return 0;
    }
  }

  // Group each color in turn.
  for (current_color = 1;
       (label_index < node_count);
       current_color++)
  {

    // Collect all nodes with current color into label.
    for (node_index = 0;
         (node_index < node_count) && (label_index < node_count);
         node_index++)
    {
      if (current_color == a->colors.a[node_index])
      {
        a->lab.a[label_index] = node_index;
        a->ptn.a[label_index++] = 1;
      }
    }

    if (label_index > 0)
    {
      // Previous label was end of current partition.
      a->ptn.a[label_index - 1] = 0;
    }
  }

  return 1;
}

/**
** Create partitions when colors are already normalized.
**
** @precondition colors are already normalized
**
** @param a object to manipulate
*/
void
cgraph_fix_ptn(cgraph_t *a)
{
  int i;
  int ip1;

  // Setup partitions assuming colors are already normalized.
  for (i = 0, ip1 = 1; ip1 < a->n; i++, ip1++)
  {
    a->lab.a[i] = i;

    if (a->colors.a[i] == a->colors.a[ip1])
    {
      // Current partition continues
      a->ptn.a[i] = 1;
    }
    else
    {
      // End of current partition
      a->ptn.a[i] = 0;
    }
  }
  // Last element
  a->lab.a[a->n - 1] = a->n - 1;
  a->ptn.a[a->n - 1] = 0;
}

/**
** Free all memory associated with cgraph.  Does not free cgraph itself.
**
** @param a object to manipulate
** @return 0 on failure, 1 otherwise
*/
boolean
cgraph_free(cgraph_t *a)
{
  int i;
  dynints_t **pcurrent_states;

  free(a->g.v);
  a->g.v = NULL;
  free(a->lab.v);
  a->lab.v = NULL;
  free(a->ptn.v);
  a->ptn.v = NULL;
  free(a->colors.v);
  a->colors.v = NULL;
  free(a->adjacent_counts.v);
  a->adjacent_counts.v = NULL;
  free(a->orbits.v);
  a->orbits.v = NULL;
  for (i = 0, pcurrent_states = a->state_dynints_array.a;
       i < a->n_alloc;
       i++, pcurrent_states++)
  {
    dynints_destroy(*pcurrent_states);
  }
  free(a->state_dynints_array.v);

  return TRUE;
}

/**
** Frees all memory associated with colored graph, including itself.
**
** @param a object to free
** @return 0 on failure, 1 otherwise
*/
boolean
cgraph_destroy(cgraph_t *a)
{
  int status;

  status = cgraph_free(a);

  free(a);

  return status;
}

/**
** Copy all data from src to new cgraph, allocating memory as needed.
**
** @param src object to copy
** @return cgraph object that is identical to src, NULL on failure
*/
cgraph_t *
cgraph_copy(cgraph_t *src)
{
  cgraph_t *dest;
  boolean status;

  dest = cgraph_create_raw(src->n);
  if (NULL == dest)
  {
    return NULL;
  }

  status = cgraph_ncopy(src, dest);
  if (FALSE == status)
  {
    cgraph_destroy(dest);
    return NULL;
  }

  return dest;
}

/**
** Copy all data from src to dest, allocating memory as needed.
**
** @param src object to copy
** @param dest object to copy to
** @return 0 on failure, 1 otherwise
*/
boolean
cgraph_ncopy(cgraph_t *src, cgraph_t *dest)
{
  int status;
  int i;
  int src_n;
  dynints_t **pcurrent_dest_states;
  dynints_t **pcurrent_src_states;

  src_n = src->n;

  status = cgraph_realloc_clear(dest, src_n);
  if (FALSE == status)
  {
    return FALSE;
  }

  memcpy(dest->g.a, src->g.a, src->m * src_n * sizeof(graph));
  memcpy(dest->lab.a, src->lab.a, src_n * sizeof(int));
  memcpy(dest->ptn.a, src->ptn.a, src_n * sizeof(int));
  memcpy(dest->colors.a, src->colors.a, src_n * sizeof(int));
  memcpy(dest->adjacent_counts.a, src->adjacent_counts.a, src_n * sizeof(int));
    
  // No need to copy orbits.

  // Copy all states.
  for (i = 0, pcurrent_dest_states = dest->state_dynints_array.a,
              pcurrent_src_states = src->state_dynints_array.a;
       i < src_n;
       i++, pcurrent_dest_states++, pcurrent_src_states++)
  {
    status = dynints_copy(*pcurrent_src_states, *pcurrent_dest_states);
    if (FALSE == status)
    {
      return FALSE;
    }
  }

  return TRUE;
}

/**
** Normalize colors and graph so cgraph can be compared with other cgraphs.
**
** @param a object to manipulate
** @return 0 on failure, 1 otherwise
*/
boolean
cgraph_normalize(cgraph_t *a)
{
  boolean status;

  // Normalize colors.
  cgraph_normalize_partitions(a);

  // Obtain canonical version of cgraph.
  status = cgraph_canon(a);

  return status;
}

/**
** Make colored graph canonical.
**
** @precondition Colors must already be normalized.
**
** @param a object to manipulate
** @return 0 on failure, 1 otherwise
*/
boolean
cgraph_canon (cgraph_t *a)
{
  static DEFAULTOPTIONS(options);
  statsblk(stats);
  set workspace[120];
  graph *canong;
  int *label_map;
  int count;
  int alloc_count;
  int i;
  dynints_t **states;
  dynints_t **old_states;
  int *colors;
  int *old_colors;
  size_t graph_alloc;

  // Must create graph with allocated size, not actual size.
  graph_alloc = ((a->n_alloc + WORDSIZE - 1) / WORDSIZE) * a->n_alloc;
  canong = calloc(graph_alloc, sizeof(graph));

  /* Generate canonically labelled isomorph for comparison */
  options.getcanon = TRUE;
  /* No directed edges or loops (default) */
  //options.digraph = FALSE;
  /* Don't output generators (default) */
  //options.writeautoms = FALSE;
  /* Don't output markers (default) */
  //options.writemarkers = FALSE;
  /* Coloring specified in ptn */
  options.defaultptn = FALSE;

  // Generate canonical graph
  nauty (a->g.a, a->lab.a, a->ptn.a, NULL, a->orbits.a, &options, &stats,
         workspace, 120, a->m, a->n, canong);

  // Get rid of old graph
  free(a->g.a);

  // Replace with new graph
  a->g.a = canong;

  // Check for errors.
  if (0 != stats.errstatus)
  {
    return FALSE;
  }

  // Use label map to reorder states and colors.
  label_map = a->lab.a;

  // Reorder states and colors to match new vertex order.
  count = a->n;
  old_states = a->state_dynints_array.a;
  old_colors = a->colors.a;
  // Must create states and colors with allocated size, not actual size.
  states = calloc((size_t)a->n_alloc, sizeof(dynints_t *));
  colors = calloc((size_t)a->n_alloc, sizeof(int));
  for (i = 0; i < count; i++)
  {
    states[i] = old_states[label_map[i]];
    colors[i] = old_colors[label_map[i]];
  }

  // Transfer remaining states.
  alloc_count = a->n_alloc;
  for (i = count; i < alloc_count; i++)
  {
    states[i] = old_states[i];
  }

  // Replace old states. Dynints are being reused, do not destroy them.
  free(old_states);
  a->state_dynints_array.a = states;

  // Replace old colors.
  free(old_colors);
  a->colors.a = colors;

  // Colors have been reordered, so recreate partitions.
  cgraph_normalize_partitions(a);

  return TRUE;
}

/**
** Compare two nodes of a graph to see if they are equivalent
** (Same color, states and adjacent to same nodes).
**
** @param a object to examine
** @param e first node to compare
** @param f second node to compare
** @return 1 when e and f are equivalent.
*/
boolean
cgraph_compare_nodes (cgraph_t *a, int e, int f)
{
  set *e_gv;
  set *f_gv;
  set *work_set;
  boolean nodes_connected;
  boolean match;
  boolean status;

  if (e == f)
  {
    // They are one and the same.
    return TRUE;
  }

  // Compare color.
  if (a->colors.a[e] != a->colors.a[f])
  {
    return FALSE;
  }

  // Compare adjacent counts.
  if (a->adjacent_counts.a[e] != a->adjacent_counts.a[f])
  {
    return FALSE;
  }

  // Compare states.
  status = dynints_compare(a->state_dynints_array.a[e],
                           a->state_dynints_array.a[f]);
  if (FALSE == status)
  {
    return FALSE;
  }

  // Compare connections.
  // If each node is connected to the same nodes, they are equivalent.
  // Ignore any connections between e and f.
  e_gv = GRAPHROW(a->g.a, e, a->m);
  f_gv = GRAPHROW(a->g.a, f, a->m);

  nodes_connected = ISELEMENT(e_gv, f);
  if (nodes_connected)
  {
    // Nodes are connected, temporarily remove connection for comparison.
    DELELEMENT(e_gv, f);
    DELELEMENT(f_gv, e);
  }

  match = (0 == memcmp(e_gv, f_gv, a->m * sizeof(set)));

  if (nodes_connected)
  {
    // Reconnect nodes.
#ifdef AGGRESSIVE_CHECK
assert(e <= a->n);
assert(f <= a->n);
#endif
    ADDELEMENT(e_gv, f);
    ADDELEMENT(f_gv, e);
  }

  if (FALSE == match)
  {
    int i;

    // Initial comparison failed, check if all connected nodes
    // are equivalent.
    work_set = malloc(a->m * sizeof(set));
    if (NULL == work_set)
    {
      error_printf("cgraph_compare_nodes can't allocate work_set, m=%d", a->m);
    }
    memcpy(work_set, f_gv, a->m * sizeof(set));

    // Don't try to match current pair.
    DELELEMENT(work_set, e);

    for (i = -1; (i = nextelement(e_gv, a->m, i)) >= 0;)
    {
      int j;

      // Don't try to match current pair.
      if (i == f)
      {
        continue;
      }

      // For each node adjacent to e, see if any unmatched node adjacent
      // to f matches.
      match = FALSE;
      for (j = -1; (j = nextelement(work_set, a->m, j)) >= 0;)
      {
        status = cgraph_compare_nodes(a, i, j);

        if (TRUE == status)
        {
          // Found match, remove from work set.
          DELELEMENT(work_set, j);

          match = TRUE;

          break;
        }
      }

      if (FALSE == status)
      {
        // No match was found for element.
        match = FALSE;
        break;
      }
    }
  }

  return match;
}

/**
** Compare two colored graphs to determine if they are equal.
**
** @precondition a and b must both be canonical and isomorphic.
**
** @param a first object to compare
** @param b second object to compare
** @return 1 when a and b are identical.
*/
boolean
cgraph_compare (cgraph_t *a, cgraph_t *b)
{
  int i;
  int a_n;
  dynints_t **pcurrent_a_states;
  dynints_t **pcurrent_b_states;

  /* Check graph size */
  /* g_sz == (m * n) */
  /* lab_sz == colors_sz == n */
  a_n = a->n;
  if ((a_n != b->n) || (a->m != b->m))
  {
    return FALSE;
  }

  /* Check graph connections */
  if (0 != memcmp(a->g.a, b->g.a, a->m * a_n * sizeof(graph)))
  {
    return FALSE;
  }

  /* Check node colors */
  if (0 != memcmp(a->colors.a, b->colors.a, a_n * sizeof(int)))
  {
    return FALSE;
  }

  /* Check node states */
  for (i = 0, pcurrent_a_states = a->state_dynints_array.a,
              pcurrent_b_states = b->state_dynints_array.a;
       i < a_n;
       i++, pcurrent_a_states++, pcurrent_b_states++)
  {
    int status;

    status = dynints_compare(*pcurrent_a_states,
                             *pcurrent_b_states);
    if (FALSE == status)
    {
      return FALSE;
    }
  }

  return TRUE;
}

/**
** Create a new colored graph c that contains everything in a and b,
** with an additional line between v_a and v_b.
**
** @param a graph to combine, must be isometric
** @param b graph to combine, must be isometric
** @param v_a vertex in graph a to join
** @param v_b vertex in graph b to join
** @param c empty graph to store the combined data, not isometric
**
** @return 0 when error occurs, 1 otherwise
*/
boolean
cgraph_link (cgraph_t *a, int v_a, cgraph_t *b, int v_b,
         cgraph_t *c)
{
  int *a_renum; // Maps element numbers from a to c
  int *b_renum; // Maps element numbers from b to c
  int status;

  a_renum = calloc((size_t)a->n, sizeof(int));
  if (NULL == a_renum) return FALSE;
  b_renum = calloc((size_t)b->n, sizeof(int));
  if (NULL == b_renum) return FALSE;

  status = cgraph_append(a, a_renum, b, b_renum, c);

  // Add link
#ifdef AGGRESSIVE_CHECK
assert(a_renum[v_a] <= c->n);
assert(b_renum[v_b] <= c->n);
#endif
  ADDELEMENT(GRAPHROW(c->g.a, b_renum[v_b], c->m), a_renum[v_a]);
  ADDELEMENT(GRAPHROW(c->g.a, a_renum[v_a], c->m), b_renum[v_b]);

  free(a_renum);
  free(b_renum);

  return TRUE;
}

/**
** Merge all vertices in two graphs into a single new graph.
**
** @param a first graph to merge
** @param a_renum array of ints of size a->n
** @param b second graph to merge
** @param b_renum array of ints of size a->n
** @param c graph that vertices from a and b will be stored
** @returns 0 on error, 1 otherwise
*/
boolean
cgraph_append (cgraph_t *a, int *a_renum,
               cgraph_t *b, int *b_renum,
               cgraph_t *c)
{
  int v;
  int i;
  int a_i;
  int b_i;
  int c_i;
  set *a_gv;
  set *b_gv;
  set *c_gv;
  boolean status;

  // Setup c.
  status = cgraph_realloc_clear(c, a->n + b->n);
  if (FALSE == status)
  {
    return FALSE;
  }

  // Copy colors and states over.
  for (a_i = 0, b_i = 0, c_i = 0; c_i < c->n; c_i++)
  {
    if (b_i == b->n)
    {
      // No more elements left in b
      a_renum[a_i] = c_i;
      c->colors.a[c_i] = a->colors.a[a_i];
      status = dynints_copy(a->state_dynints_array.a[a_i],
                            c->state_dynints_array.a[c_i]);
      if (FALSE == status)
      {
        return FALSE;
      }
      a_i++;
    }
    else if (a_i == a->n)
    {
      // No more elements left in a
      b_renum[b_i] = c_i;
      c->colors.a[c_i] = b->colors.a[b_i];
      status = dynints_copy(b->state_dynints_array.a[b_i],
                            c->state_dynints_array.a[c_i]);
      if (FALSE == status)
      {
        return FALSE;
      }
      b_i++;
    }
    else if (a->colors.a[a_i] <= b->colors.a[b_i])
    {
      // Need to maintain colors in sorted order.
      a_renum[a_i] = c_i;
      c->colors.a[c_i] = a->colors.a[a_i];
      status = dynints_copy(a->state_dynints_array.a[a_i],
                            c->state_dynints_array.a[c_i]);
      if (FALSE == status)
      {
        return FALSE;
      }
      a_i++;
    }
    else
    {
      b_renum[b_i] = c_i;
      c->colors.a[c_i] = b->colors.a[b_i];
      status = dynints_copy(b->state_dynints_array.a[b_i],
                            c->state_dynints_array.a[c_i]);
      if (FALSE == status)
      {
        return FALSE;
      }
      b_i++;
    }
  }

  // Insert a's elements into c
  for (v = 0;
       v < a->n;
       v++)
  {
    a_gv = GRAPHROW(a->g.a, v, a->m);
    c_gv = GRAPHROW(c->g.a, a_renum[v], c->m);
    for (i = -1; (i = nextelement(a_gv, a->m, i)) >= 0;)
    {
#ifdef AGGRESSIVE_CHECK
assert(a_renum[i] <= c->n);
#endif
      ADDELEMENT(c_gv, a_renum[i]);
    }
  }

  // Insert b's elements into c
  for (v = 0;
       v < b->n;
       v++)
  {
    b_gv = GRAPHROW(b->g.a, v, b->m);
    c_gv = GRAPHROW(c->g.a, b_renum[v], c->m);
    for (i = -1; (i = nextelement(b_gv, b->m, i)) >= 0;)
    {
#ifdef AGGRESSIVE_CHECK
assert(b_renum[i] <= c->n);
#endif
      ADDELEMENT(c_gv, b_renum[i]);
    }
  }

  // Recreate partitions.
  cgraph_fix_ptn(c);

  return TRUE;
}

/**
** Create two new colored graphs c and d that contain
** the seperate graphs created by deleting the link between
** v_a and v_b.
**
** NOTE: d->n will be 0 if removing link does not create two separate
**       graphs, i.e. all vertices in graph a are still connected in some way.
**
** @postcondition c is modified to contain the subgraph containing v_a.
** @postcondition d is modified to contain the subgraph containing v_b.
** @postcondition c and d are NOT isometric after calling this function.
**
** @param a graph to delete line from, must be isometric.
** @param v_a vertex in graph a to unlink.
** @param v_b vertex in graph b to unlink.
** @param c empty graph to store the subgraph containing v_a.
** @param d empty graph to store the subgraph containing v_b.
**
** @return 0 on error, 1 otherwise
*/
boolean
cgraph_unlink (cgraph_t *a, int v_a, int v_b,
               cgraph_t *c, cgraph_t *d)
{
  boolean status;

//putgraph(custom_get_fh(ptype), a->g.a, 80, a->m, a->n);
//i=0;
//putset(custom_get_fh(ptype), GRAPHROW(a->g.a, v_a, a->m), &i, 80, a->m, FALSE);
//printf("\n");
  // Break link between v_a and v_b.
  DELELEMENT(GRAPHROW(a->g.a, v_a, a->m), v_b);
  DELELEMENT(GRAPHROW(a->g.a, v_b, a->m), v_a);

  status = cgraph_seperate(a, c, d);

  // Restore links that were deleted from a.
#ifdef AGGRESSIVE_CHECK
assert(v_a <= a->n);
assert(v_b <= a->n);
#endif
  ADDELEMENT(GRAPHROW(a->g.a, v_a, a->m), v_b);
  ADDELEMENT(GRAPHROW(a->g.a, v_b, a->m), v_a);

  return status;
}

/**
** Seperate graph if it is not fully connected.
** 
** NOTE: d->n will be 0 if removing link does not create two separate
**       graphs, i.e. all vertices in graph a are still connected in some way.
** NOTE: If c is NULL, no graph operations are performed
**
** @param a input graph
** @param c output graph, contains fully connected graph
** @param d output graph, contains remaining graph(s), if any
** @return 0 on failure, 1 otherwise; if c is NULL, 2 if split, 0 on failure
**
** @postcondition c is modified to contain the subgraph containing v_a
** @postcondition d is modified to contain the subgraph containing v_b
** @postcondition c and d are NOT isometric after calling this function
*/
boolean
cgraph_seperate (cgraph_t *a, cgraph_t *c, cgraph_t *d)
{
  boolean status;
  int *c_renum; // Maps element numbers of A to C.
  int *d_renum; // Maps element numbers of A to D.

  // c->m and d->m need to be minimal to facilitate matching.

  c_renum = calloc((size_t)a->n, sizeof(int));
  d_renum = calloc((size_t)a->n, sizeof(int));

//putset(custom_get_fh(ptype), GRAPHROW(a->g.a, v_a, a->m), &i, 80, a->m, FALSE);
//printf("\n");

//putgraph(custom_get_fh(ptype), a->g.a, 80, a->m, a->n);

  status = cgraph_seperate_raw(a, c, d, c_renum, d_renum);

  free(c_renum);
  free(d_renum);

  return status;
}

/**
** Seperate graph if it is not fully connected.
** 
** NOTE: d->n will be 0 if removing link does not create two separate
**       graphs, i.e. all vertices in graph a are still connected in some way.
** NOTE: If c is NULL, no graph operations are performed
**
** @param a input graph
** @param c output graph, contains fully connected graph
** @param d output graph, contains remaining graph(s), if any
** @param c_renum array of ints of size a->n
** @param d_renum array of ints of size a->n
** @return 0 on failure, 1 otherwise; if c is NULL, 2 if split, 1 otherwise
**
** @postcondition c is modified to contain the subgraph containing v_a.
** @postcondition d is modified to contain the subgraph containing v_b.
** @postcondition c and d are NOT isometric after calling this function.
*/
boolean
cgraph_seperate_raw (cgraph_t *a, cgraph_t *c, cgraph_t *d,
                     int *c_renum, int *d_renum)
{
  boolean status;
  set *work_set; // List of vertices remaining to explore.
  set *c_set; // List of all vertices in subgraph containing v_a.

  work_set = calloc((size_t)a->m, sizeof(set));
  c_set = calloc((size_t)a->m, sizeof(set));

  status = cgraph_seperate_raw_raw(a, c, d, work_set, c_set, c_renum, d_renum);

  free(work_set);
  free(c_set);

  return status;
}

/**
** Seperate graph if it is not fully connected.
** 
** NOTE: d->n will be 0 if removing link does not create two separate
**       graphs, i.e. all vertices in graph a are still connected in some way.
** NOTE: If c is NULL, no graph operations are performed
**
**
** @param a input graph
** @param c output graph, contains fully connected graph
** @param d output graph, contains remaining graph(s), if any
** @param work_set set object of size m.
** @param c_set set object of size m.
** @param c_renum array of ints of size a->n
** @param d_renum array of ints of size a->n, set to NULL if graphs left alone.
** @return 0 on failure, 1 otherwise; if c is NULL, 2 if split, 1 otherwise
**
** @postcondition c is modified to contain the subgraph containing v_a.
** @postcondition d is modified to contain the subgraph containing v_b.
** @postcondition c and d are NOT isometric after calling this function.
*/
static boolean
cgraph_seperate_raw_raw (cgraph_t *a, cgraph_t *c, cgraph_t *d,
                         set *work_set, set *c_set, int *c_renum, int *d_renum)
{
  set *work2_set; // List of vertices current vertex is connected to.
  int next;
  int i;
  int j;
  int v;
  set *a_gv;
  set *c_gv;
  set *d_gv;
  boolean status;
  int first_element;

  first_element = 0;

  // Initialize c_set to contain first_element;
  EMPTYSET(c_set, a->m);
  ADDELEMENT(c_set, first_element);

  // Initialize work_set to contain first_element.
  EMPTYSET(work_set, a->m);
  ADDELEMENT(work_set, first_element);

  // Collect vertex list for graph C.

  // Loop through each element until no more elements to check.
  for (next = -1; (next = nextelement(work_set, a->m, next)) >= 0;)
  {
    int next2;

    // Remove current element from work_set so it is not checked again.
    DELELEMENT(work_set, next);

    // Get list of elements current element is connected to.
    work2_set = GRAPHROW(a->g.a, next, a->m);

    // Loop through each element connected to current element
    for (next2 = -1; (next2 = nextelement(work2_set, a->m, next2)) >= 0;)
    {
      if (!ISELEMENT(c_set, next2))
      {
        // Not part of c_set yet, so add to work_set and c_set.
#ifdef AGGRESSIVE_CHECK
assert(next2 <= a->n);
#endif
        ADDELEMENT(work_set, next2);
        ADDELEMENT(c_set, next2);

        // Restart outer loop because new elements added
        next = -1;
      }
    }
  }

  // Specify renumbering for all elements in c.
  i = 0;
  j = 0;
  for (next = 0; next < a->n; next++)
  {
    if (ISELEMENT(c_set, next))
    {
      // Connected to first_element.
      c_renum[next] = i;
      d_renum[next] = -1;
      i++;
    }
    else
    {
      // Not connected to first_element.
      d_renum[next] = j;
      c_renum[next] = -1;
      j++;
    }
  }

  if (NULL == c)
  {
    // Don't modify graph, simply exit.
    if (0 == j)
    {
      // Graph was not split.
      return 1;
    }
    else
    {
      // Graph was split.
      return 2;
    }
  }

  if (0 == j)
  {
    // Graph did not get split, so we can simply copy graph with missing link.
    status = cgraph_ncopy(a, c);
    if (FALSE == status)
    {
      return FALSE;
    }

    // Mark d as empty.
    status = cgraph_realloc_clear(d, 0);
    if (FALSE == status)
    {
      return FALSE;
    }
  }
  else
  {
    // Graph got split in two.

    //printf("i = %i, j = %i\n", i, j);
    // Setup c.
    status = cgraph_realloc_clear(c, i);
    if (FALSE == status)
    {
      return FALSE;
    }

    // Setup d.
    status = cgraph_realloc_clear(d, j);
    if (FALSE == status)
    {
      return FALSE;
    }


    // Copy elements to new graphs.
    for (v = 0, a_gv = a->g.a; v < a->n; v++, a_gv += a->m)
    {
      if (ISELEMENT(c_set, v))
      {
        // Put element into graph c.

        // Get appropriate element set
        c_gv = GRAPHROW(c->g.a, c_renum[v], c->m);

        // Recreate connections using renumbered elements
        for (next = -1; (next = nextelement(a_gv, a->m, next)) >= 0;)
        {
#ifdef AGGRESSIVE_CHECK
assert(c_renum[next] <= c->n);
#endif
          ADDELEMENT(c_gv, c_renum[next]);
        }
      }
      else
      {
        // Put element into graph d.

        // Get appropriate element set
        d_gv = GRAPHROW(d->g.a, d_renum[v], d->m);

        // Recreate connections using renumbered elements
        for (next = -1; (next = nextelement(a_gv, a->m, next)) >= 0;)
        {
#ifdef AGGRESSIVE_CHECK
assert(d_renum[next] <= d->n);
#endif
          ADDELEMENT(d_gv, d_renum[next]);
        }
      }
    }

    // Copy color of each element to new graphs.
    for (v = 0; v < a->n; v++)
    {
      if (ISELEMENT(c_set, v))
      {
        // Element in graph c.
        c->colors.a[c_renum[v]] = a->colors.a[v];
        status = dynints_copy(a->state_dynints_array.a[v],
                              c->state_dynints_array.a[c_renum[v]]);
        if (FALSE == status)
        {
          return FALSE;
        }
        c->lab.a[c_renum[v]] = a->lab.a[v];
      }
      else
      {
        // Element in graph d.
        d->colors.a[d_renum[v]] = a->colors.a[v];
        status = dynints_copy(a->state_dynints_array.a[v],
                              d->state_dynints_array.a[d_renum[v]]);
        if (FALSE == status)
        {
          return FALSE;
        }
        d->lab.a[d_renum[v]] = a->lab.a[v];
      }
    }

    // Create partitions for c
    cgraph_fix_ptn(c);

    // Create partitions for d
    cgraph_fix_ptn(d);

  }

  return TRUE;
}

