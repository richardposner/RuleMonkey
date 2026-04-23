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


#ifndef PARSER_LIB_H
#define PARSER_LIB_H

typedef struct bng_command_args_struct bng_command_args_t;
typedef enum block_type_enum block_type_t;

#include "constants.h"
#include "dynsymbols.h"
#include "dyncolors.h"
#include "dynmolecules.h"
#include "dynconnections.h"
#include "dynmaps.h"
#include "dynspecies.h"
#include "cgraph.h"
#include "reactant_definition.h"

struct bng_command_args_struct
{
  // Base name of CDAT and GDAT files.
  // Default: Base name of BNGL file.
  char *prefix;


  // Absolute error tolerance for species concentrations.
  // Default: 1e-8.
  double atol;

  // Relative error tolerance for species concentrations.
  // Default: 1e-8.
  double rtol;

  // Turns on use of sparse matrix formation of the Jabian and
  // iterative solution of linear equations using GNRES.
  // Recommended for networks with more than a few hundred species.
  // Default: 0 (Off).
  int sparse;


  // Number of intervals at which to report concentrations.
  // Default: 1.
  double n_steps;

  // Starting time for integration.
  // Default: 0
  double t_start;

  // End time for integration.
  // Default: None.
  double t_end;

  // Seed for random number generator.
  // Default: 0.
  unsigned seed;

  // Volume of space that is involved in the simulation.
  // Default: 1.
  double volume;

  // Size of timeslice, used for particle based simulations.
  // Default: 0 (will be computed automatically if not specified).
  double time_delta;

  // Maximum value that a computed probability is allowed to have before
  // time_delta is automatically recomputed.
  // Default: 0.9
  double max_prob;

  // When != 0, multiply the number of particles by pseudo_factor
  // to get the number of pseudo particles.
  double pseudo_factor;

  // When != 0, forces number of pseudo particles to pseudo_count.
  double pseudo_count;

  // String to add to the output filenames.
  char *suffix;
};

enum block_type_enum
{
  NO_BLOCK_TYPE,
  PARAMETERS_BLOCK_TYPE,
  SPECIES_BLOCK_TYPE,
  REACTION_RULES_BLOCK_TYPE,
  OBSERVABLES_BLOCK_TYPE,
  MOLECULE_TYPES_BLOCK_TYPE,
  ACTIONS_BLOCK_TYPE
};

// Global variables.
/* Contains the simulated world according to dynstoc. */
extern world_t *world;

/* If TRUE, enforce stricter syntax */
extern int strict;

/* Symbol lookup table. */
extern dynsymbols_t *dynsymbols;

/* Molecules, used to validate configuration file. */
extern dynmolecules_t *dynmolecules;

/* Associates molecules on lhs of reaction with molecules on rhs. */
extern dynmaps_t *dynmaps;

/* Current molecule being operated on. */
extern molecule_t *current_molecule;
extern int current_molecule_index;

/* Current molecule component being operated on. */
extern int current_component_index;

/* Object index, refers to molecule or component depending on context. */
extern int current_object_index;
extern int current_object_rule_index;
extern const char *current_object_name;

/* Used to build rules. */
extern int current_molecule_rule_index;

/* Used to verify current molecule. */
extern int *molecule_component_multipliers_table;

/* Used to build specie definition during specie initialization. */
extern cgraph_t *current_specie_cgraph;

/* Used to keep track of connections between molecules. */
extern dynconnections_t *dynconnections;

/* TRUE only if molecule didn't exist before. */
extern boolean new_molecule;

/* TRUE only if initializing new specie (or resetting population) */
extern boolean initializing_specie;

/* TRUE only if building specie matching rule. */
extern boolean building_specie;

/* TRUE only if building new rule */
extern boolean building_rule;

/* Number of reaction definitions currently defined.  */
extern int reaction_definition_count;

/* For more descriptive error messages. */
extern const char *current_molecule_name;

/* To keep track of all the reactant/product rule definitions. */
extern boolean reaction_lhs;
extern int specie_rule_indexes[MAX_SPECIE_RULE_INDEX];
extern int current_specie_index;

/* Information needed to run BNG command. */
extern boolean (*bng_command_ptr)(bng_command_args_t *);
extern bng_command_args_t bng_command_args;

/* So that name of reaction definition can be set. */
extern reaction_definition_t *forward_reaction_definition;
extern reaction_definition_t *reverse_reaction_definition;

/* Keep track of what block we are currently processing. */
extern block_type_t current_block;

/* Reinforce molecule type definitions if present. */
extern int molecule_types_present;

/* TRUE if '*' is present in reaction/product definition. */
extern int molecule_star_count;
extern int molecule_star_lhs;
extern int molecule_star_rhs;

// yyerror defined in ds_bngl.l.
void yyerror(char *s);
void yyerror_prev(char *s);

boolean setup_parse_env();
int setup_molecule_name(const char *color_name, const boolean initialize);
boolean add_object_modifier(const color_type_t type, const char *name);
boolean add_connection(const int specie_rule_index,
                       const char *name, const int node);
boolean finish_reaction_rules(const char *name, boolean bidirectional,
                              const int lhs_count, const int rhs_count,
                              double forward_rate, double reverse_rate);
boolean finish_observable(const char *name, const int rule_index,
                          double denominator);
boolean setup_reaction_rule(const char *name, const double rate,
                            const int flags, const int reactant_count,
                            const int product_count,
                            int *current_specie_rule_indexes,
                            dynints_t *new_species);
boolean bng_command_set_default();
boolean bng_command_set_double(const char *name, const double value);
boolean bng_command_set_string(const char *name, const char *value);
boolean setup_ps(bng_command_args_t *args);
boolean simulate_ps(bng_command_args_t *args);
boolean simulate_rs(bng_command_args_t *args);
boolean state_output_ps(bng_command_args_t *args);
boolean ds_system_set_value(const char *name, double value);
boolean ds_system_get_value(const char *name, double *pvalue);


#endif /* DS_BNGL_H */
