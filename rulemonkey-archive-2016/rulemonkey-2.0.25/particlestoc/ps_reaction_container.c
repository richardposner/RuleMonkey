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


#include <assert.h>
#include "ps_reaction_container.h"
#include "dynps_species.h"
#include "rule_aggregate.h"

void suppress_warnings26() {SUPPRESS_WARNINGS;}

void
ps_reaction_container_setup(const size_t full_aggregate_size,
                            reaction_definition_t *reaction_definition,
                            const int reactant_index1, int *aggregate1,
                            const int reactant_index2, int *aggregate2,
                            ps_reaction_container_t *ps_reaction_container)
{
  ps_reaction_container->reaction_definition = reaction_definition;
  ps_reaction_container->reactant_index1 = reactant_index1;
  ps_reaction_container->reactant_index2 = reactant_index2;
  ps_reaction_container->aggregate1 = aggregate1;
  ps_reaction_container->aggregate2 = aggregate2;
  ps_reaction_container->full_aggregate_size = full_aggregate_size;
}

boolean
ps_reaction_container_perform_reaction
  (ps_reaction_container_t *ps_reaction_container,
   ps_particle_t *particle1, ps_particle_t *particle2,
   world_t *world, dynps_species_t *dynps_species)
{
  ps_particle_t *full_particle;
  cgraph_t *final_cgraph;
  int *full_aggregate;
  int full_aggregate_allocated;
  specie_t *specie1;
  specie_t *specie2;
  int status;

  // Count number of reactions done.
  world->reaction_count += 1;

  specie1 = world->dynspecies->species.a + particle1->specie_index;
  specie2 = NULL;

  if ((NULL == particle2) ||
      (ps_particle_check_same_entity(particle1, particle2)))
  {
    // Unimolecular (uniparticle or biparticle) reaction, nothing to merge.
    full_particle = particle1;
    final_cgraph = cgraph_copy(specie1->cg);
    full_aggregate = ps_reaction_container->aggregate1;
    full_aggregate_allocated = 0;

    // Update population, but don't remove specie until after reaction has
    // occured so that we can use aggregate.
    specie_update(specie1, -1, world);
  }
  else
  {
    // Biparticle reaction.
    specie2 = world->dynspecies->species.a + particle2->specie_index;
  
    // Will be initialized by ps_particle_cgraph_append.
    full_particle = NULL;

    // Create cgraph to store combined cgraph.
    final_cgraph = cgraph_create(INITIAL_CGRAPH_SIZE);

    // Create memory to merge aggregates into.
    full_aggregate = calloc(sizeof(int),
                            ps_reaction_container->full_aggregate_size);
    if (NULL == full_aggregate)
    {
      assert(0);
      return FALSE;
    }
    full_aggregate_allocated = 1;

    // Add delimiters to initialize full aggregate.
    rule_aggregate_init_full
      (full_aggregate, ps_reaction_container->aggregate1);

    // Merge first aggregate into full_aggregate.
    // Also copy particle1 to full_particle.
    status = ps_particle_cgraph_append(&full_particle, NULL,
                                       full_aggregate,
                                       particle1, NULL,
                                       ps_reaction_container->aggregate1,
                                       NULL);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "Merged 1 to full_aggregate:\n");
    RULE_AGGREGATE_OUTPUT(OUTPUT_DEBUG, 5, ps_reaction_container->aggregate1);
    DEBUG_PRINTF(5, "\n");
    RULE_AGGREGATE_OUTPUT(OUTPUT_DEBUG, 5, full_aggregate);
    DEBUG_PRINTF(5, "\n");
#endif

    // Merge second aggregate into full aggregate and merge both cgraphs.
    // Also merge particle1 into particle2, renumbering appropriately.
    status = ps_particle_cgraph_append(&full_particle, specie1->cg,
                                       full_aggregate,
                                       particle2, specie2->cg,
                                       ps_reaction_container->aggregate2,
                                       final_cgraph);
    if (FALSE == status)
    {
      assert(0);
      return FALSE;
    }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "Merged 2 to full_aggregate:\n");
    RULE_AGGREGATE_OUTPUT(OUTPUT_DEBUG, 5, ps_reaction_container->aggregate2);
    DEBUG_PRINTF(5, "\n");
    RULE_AGGREGATE_OUTPUT(OUTPUT_DEBUG, 5, full_aggregate);
    DEBUG_PRINTF(5, "\n");
#endif

    // Update populations, but don't remove species until after reaction has
    // occured so that we can use aggregate.
    specie_update(specie1, -1, world);
    specie_update(specie2, -1, world);
  }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "Before applying aggregate:\n");
    CGRAPH_OUTPUT(OUTPUT_DEBUG, 5, final_cgraph, world->dyncolors);
#endif

  // Perform reaction.
  status = rule_aggregate_cgraph_apply(final_cgraph, full_aggregate);
  if (FALSE == status)
  {
    assert(0);
    if (full_aggregate_allocated != 0)
    {
      free(full_aggregate);
    }
    return FALSE;
  }

  // We can remove species if needed now that reaction has finished.
  if ((specie2 != NULL) && (0 == specie2->pop) && (specie1->id != specie2->id))
  {
    dynps_species_specie_remove(specie2, dynps_species);
  }
  if (0 == specie1->pop)
  {
    dynps_species_specie_remove(specie1, dynps_species);
  }

#ifdef RULE_CHECK_DEBUG
    DEBUG_PRINTF(5, "After applying aggregate:\n");
    CGRAPH_OUTPUT(OUTPUT_DEBUG, 5, final_cgraph, world->dyncolors);
#endif

  // Add product(s) to population, splitting cgraph if needed.
  status = ps_particle_add_new_species(full_particle, final_cgraph, world,
                                       dynps_species);
  if (FALSE == status)
  {
    if (full_aggregate_allocated != 0)
    {
      free(full_aggregate);
    }
    return FALSE;
  }

  cgraph_destroy(final_cgraph);

  if (full_aggregate_allocated != 0)
  {
    free(full_aggregate);
  }

  return TRUE;
}

void
ps_reaction_container_output(const output_type ptype, const int level,
                             ps_reaction_container_t *ps_reaction_container)
{
  CUSTOM_PRINTF(ptype, level, "definition: %s [",
                ps_reaction_container->reaction_definition->name);

  if (NULL != ps_reaction_container->aggregate1)
  {
    rule_aggregate_output(ptype, level, ps_reaction_container->aggregate1);
  }
  else
  {
    CUSTOM_PRINTF(ptype, level, "NULL");
  }

  CUSTOM_PRINTF(ptype, level, "] [");

  if (NULL != ps_reaction_container->aggregate2)
  {
    rule_aggregate_output(ptype, level, ps_reaction_container->aggregate2);
  }
  else
  {
    CUSTOM_PRINTF(ptype, level, "NULL");
  }

  CUSTOM_PRINTF(ptype, level, "]\n");
}
