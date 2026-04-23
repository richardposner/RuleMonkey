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
#include "particlestoc.h"
#include "myrand.h"
#include "cgraph.h"
#include "output.h"
#include <assert.h>
#include <stdio.h>
#include "reactionstoc.h"

#define FORWARD_SUFFIX "_f"
#define REVERSE_SUFFIX "_r"

#define NAME_STR_LENGTH 512

void suppress_warnings31() {SUPPRESS_WARNINGS;}

// Global variables.
/* Contains the simulated world according to rulemonkey. */
world_t *world;

// Various system variables that can be manipulated by script directly.
struct ds_system_struct
{
  int strict;

  // Constant to use for the avogadro number
  double avogadro;
};
typedef struct ds_system_struct ds_system_t;

ds_system_t ds_system;

/* Symbol lookup table. */
dynsymbols_t *dynsymbols;

/* Molecules, used to validate configuration file. */
dynmolecules_t *dynmolecules;

/* Associates molecules on lhs of reaction with molecules on rhs. */
dynmaps_t *dynmaps;

/* Current molecule being operated on. */
molecule_t *current_molecule;
int current_molecule_index;

/* Current molecule component being operated on. */
int current_component_index;

/* Object index, refers to molecule or component depending on context. */
int current_object_index;
int current_object_rule_index;
const char *current_object_name;

/* Used to build rules. */
int current_molecule_rule_index;

/* Used to verify current molecule. */
int *molecule_component_multipliers_table;

/* Used to build specie definition during specie initialization. */
cgraph_t *current_specie_cgraph;

/* Used to keep track of connections between molecules. */
dynconnections_t *dynconnections;

/* TRUE only if molecule didn't exist before. */
boolean new_molecule;

/* TRUE only if initializing new specie (or resetting population) */
boolean initializing_specie;

/* TRUE only if building specie matching rule. */
boolean building_specie;

/* TRUE only if building new rule. */
boolean building_rule;

/* Number of reaction definitions currently defined.  */
int reaction_definition_count;

/* For more descriptive error messages. */
const char *current_molecule_name;

/* To keep track of all the reactant/product rule definitions. */
boolean reaction_lhs;
int specie_rule_indexes[MAX_SPECIE_RULE_INDEX];
int current_specie_index;

/* Keep track of current block we are processing. */
block_type_t current_block;

/* Reinforce molecule definitions if molecule types present. */
int molecule_types_present;

/* Information needed to run BNG command. */
boolean (*bng_command_ptr)(bng_command_args_t *);
bng_command_args_t bng_command_args;

/* TRUE if '*' is present in reaction/product definition. */
int molecule_star_count;
int molecule_star_lhs;
int molecule_star_rhs;

static void ds_system_setup();

boolean
setup_parse_env()
{
  world = world_create();
  if (NULL == world)
  {
    yyerror("Error creating rulemonkey world");
    return FALSE;
  }

  ds_system_setup();

  dynsymbols = dynsymbols_create(INITIAL_DYNSYMBOLS);
  if (NULL == dynsymbols)
  {
    yyerror("Error creating dynsymbols");
    return FALSE;
  }

  /* For verifying molecules match any previous definitions. */
  dynmolecules = dynmolecules_create(INITIAL_DYNMOLECULES);
  if (NULL == dynmolecules)
  {
    yyerror("Error creating dynmolecules");
    return FALSE;
  }

  /* For associating molecules between lhs and rhs of reaction. */
  dynmaps = dynmaps_create(INITIAL_DYNMAPS);
  if (NULL == dynmaps)
  { 
    yyerror("Unable to create dynmaps");
    exit(1);
  }
    
  
  /* No molecules currently being operated on. */
  current_molecule = NULL;
  current_molecule_index = -1;
  current_component_index = -1;
  current_object_index = -1;
  current_object_rule_index = -1;
  current_object_name = NULL;
  current_molecule_rule_index = -1;
  molecule_component_multipliers_table = NULL;

  /* For building specie currently being initialized. */
  current_specie_cgraph = cgraph_create(INITIAL_CGRAPH_SIZE);
  if (NULL == current_specie_cgraph)
  {
    yyerror("Unable to create current_specie");
    exit(1);
  }

  /* For keeping track of inter-molecule connections. */
  dynconnections = dynconnections_create(INITIAL_DYNCONNECTIONS);
  if (NULL == dynconnections)
  {
    yyerror("Unable to create dynconnections");
    exit(1);
  }

  new_molecule = FALSE;
  initializing_specie = FALSE;
  building_specie = FALSE;
  building_rule = FALSE;
  reaction_definition_count = 0;
  current_molecule_name = NULL;
  reaction_lhs = FALSE;
  memset(specie_rule_indexes, 0, sizeof(int) * MAX_SPECIE_RULE_INDEX);
  current_specie_index = -1;
  bng_command_ptr = NULL;
  memset(&bng_command_args, 0, sizeof(bng_command_args));

  return TRUE;
} 


int
setup_molecule_name(const char *color_name, const boolean initialize)
{   
  color_t *molecule_color;

  /* Get color associated with molecule. */
  molecule_color = dyncolors_color_create(color_name, COLOR_TYPE_MOLECULE,
                                          world->dyncolors);
  if (NULL == molecule_color)
  {
    yyerror("Unable to create molecule color"); 
    exit(1);
  }

  /* Setup current_molecule. */
  current_molecule = dynmolecules_molecule_find(molecule_color->value,
                                                dynmolecules);
  if (NULL != current_molecule)
  {
    /*
     * Create temporary molecule component multipliers table
     * for comparison purposes.
     */
    assert(NULL == molecule_component_multipliers_table);
    molecule_component_multipliers_table =
      molecule_component_multipliers_table_create(current_molecule);

    new_molecule = FALSE;
  }
  else
  {
    if (FALSE == initialize)
    {
      /* Molecule must be initialized before being used. */
      yyerror("Molecule used before getting initialized");
      exit(1);
    }

    /* Create new molecule. */
    current_molecule =
      dynmolecules_molecule_create(molecule_color->value,
                                   INITIAL_DYNMOLECULES,
                                   dynmolecules);
    if (NULL == current_molecule)
    {
      yyerror("Unable to create molecule");
      exit(1);
    }

    new_molecule = TRUE;
  }

  return molecule_color->value;
}

boolean
add_object_modifier(const color_type_t type, const char *name)
{
  color_t *modifier_color;

  /* Lookup modifier color. */
  modifier_color = real_dyncolors_color_create(name, type, world->dynspecies,
                                               world->dyncolors, world);
  if (NULL == modifier_color)
  {
    yyerror("Error initializing modifier");
    return FALSE;
  }

  if (TRUE == initializing_specie)
  {
    int modifier_index;

    /* Add modifier node. */
    modifier_index = cgraph_add_node(current_specie_cgraph,
                                     modifier_color->value);
    if (-1 == modifier_index)
    {
      yyerror("Unable to add modifier to component");
      return FALSE;
    }
    /* Connect modifier to component. */
    cgraph_connect_nodes(current_specie_cgraph, current_object_index,
                         modifier_index);
  }
  else
  {
    int status;
    int modifier_rule_index;

    /* Create modifier rule into object_rule's adjacent chain. */
    status = dynrules_rule_create(modifier_color->value,
                                  ADJACENT_EXCLUSIVE_MODIFIER,
                                  EQUAL, 1, 1.0, -1, -1, -1, -1,
                                  world->dynrules);
    if (FALSE == status)
    {
      yyerror("Unable to create rule for molecule");
      return FALSE;
    }

    /* Get index of new rule. */
    modifier_rule_index = world->dynrules->count - 1;

    /* Save pointers for convenience. */
    dynrules_rule_add_adjacent_rule
              (current_object_rule_index,modifier_rule_index,
               world->dynrules);
  }

  return TRUE;
}

boolean
add_connection(const int specie_rule_index, const char *name, const int node)
{
  boolean status;

  status = dynconnections_connection_add(specie_rule_index, name, node,
                                         dynconnections);

  return status;
}

boolean
finish_reaction_rules(const char *name, boolean bidirectional,
                      const int in_lhs_count, const int in_rhs_count,
                      double forward_rate, double reverse_rate)
{
  int lhs_count;
  int rhs_count;
  char *forward_name;
  char *reverse_name;
  int flags;
  int status;
  int name_len;
  dynints_t *extra_lhs_molecules;
  dynints_t *extra_rhs_molecules;
  dynints_t *extra_lhs_species;
  dynints_t *extra_rhs_species;

  lhs_count = in_lhs_count;
  rhs_count = in_rhs_count;

  // Destroyed before function exits.
  extra_lhs_molecules = dynints_create(1);
  extra_rhs_molecules = dynints_create(1);

  // If populated, is part of new reaction rule.
  extra_lhs_species = dynints_create(1);
  extra_rhs_species = dynints_create(1);

  /* The flags may be used at some point. */
  flags = 0;

#ifdef RULE_CREATE_DEBUG
  DEBUG_PRINTF(1,
               "reaction rule: bi = %d lhs = %d, rhs = %d; rates = %f, %f\n",
               bidirectional, lhs_count, rhs_count, forward_rate,
               reverse_rate);
#endif

  // Increment reaction definition counter.
  reaction_definition_count++;

#ifdef RULE_CREATE_DEBUG
  DYNRULES_OUTPUT(OUTPUT_DEBUG, 10, world->dynrules, world->dyncolors);
  DYNCONNECTIONS_OUTPUT(OUTPUT_DEBUG, 3, dynconnections);
#endif

  if ((0 != lhs_count) && (0 != rhs_count))
  {
    // Map molecules on lhs with molecules on rhs.
    status = dynrules_map_molecules(lhs_count, rhs_count,
                                    molecule_star_lhs, molecule_star_rhs,
                                    specie_rule_indexes, dynmaps,
                                    extra_lhs_molecules, extra_rhs_molecules,
                                    dynconnections, world->dyncolors,
                                    world->dynrules);
    if (FALSE == status)
    {
      yyerror_prev("Error mapping molecules across reaction");
      dynints_destroy(extra_lhs_molecules);
      dynints_destroy(extra_rhs_molecules);
      return FALSE;
    }

#ifdef RULE_CREATE_DEBUG
    DYNMAPS_OUTPUT(OUTPUT_DEBUG, 3, dynmaps);
#endif

    // Map any split directives to corresponding binds.
    status =
      dynconnections_apply_map(lhs_count, rhs_count, specie_rule_indexes,
                               world->dynrules, dynmaps, dynconnections);
    if (FALSE == status)
    {
      yyerror_prev("Error applying map  across reaction");
      dynints_destroy(extra_lhs_molecules);
      dynints_destroy(extra_rhs_molecules);
      return FALSE;
    }

    // Use maps to discover any changed states and add directives accordingly.
    status = dynmaps_apply_to_rules(dynmaps, world->dynrules);
    if (FALSE == status)
    {
      yyerror_prev("Error applying maps to rules");
      dynints_destroy(extra_lhs_molecules);
      dynints_destroy(extra_rhs_molecules);
      assert(0);
      return FALSE;
    }

    // Apply connections between molecules and add any split directives.
    // Modifies rule hierarchy, so must be performed after map is found.
    status =
      dynconnections_apply_to_rules(lhs_count, rhs_count, ds_system.strict,
                                    specie_rule_indexes, world->dyncolors,
                                    world->dynrules, dynconnections);
    if (FALSE == status)
    {
      yyerror_prev("Error applying connections to rules");
      dynints_destroy(extra_lhs_molecules);
      dynints_destroy(extra_rhs_molecules);
      return FALSE;
    }
  }
  else
  {
    if (lhs_count != 0)
    {
      yyerror_prev("Not supported yet");
      return FALSE;
    }
    else if (rhs_count != 0)
    {
      yyerror_prev("Not supported yet");
      return FALSE;
    }
    else
    {
      yyerror_prev("No reactants or products");
      return FALSE;
    }
  }


  // Current specie creation and deletion is very limited.
  if ((0 != extra_lhs_molecules->count) || (extra_rhs_molecules->count > 1))
  {
    error_printf("extra_lhs_molecules: %zu, extra_rhs_molecules: %zu\n",
                 extra_lhs_molecules->count, extra_rhs_molecules->count);
    yyerror_prev("Error, unexpected extra molecules");
    return FALSE;
  }

  /* Modify name to distinguish between forward and reverse. */
  if ((NULL == name) || ('\0' == *name))
  {
    // Create new name based on the reaction definition number.
#ifdef HAVE_ASPRINTF
    asprintf(&forward_name, "%d_f", reaction_definition_count);
    if (NULL == forward_name) return FALSE;
    asprintf(&reverse_name, "%d_r", reaction_definition_count);
    if (NULL == reverse_name) return FALSE;
#else
    forward_name = malloc(NAME_STR_LENGTH);
    if (NULL == forward_name) return FALSE;
    reverse_name = malloc(NAME_STR_LENGTH);
    if (NULL == reverse_name) return FALSE;
    status = sprintf(forward_name, "%d_f", reaction_definition_count);
    if (status < 0) return FALSE;
    if (status > NAME_STR_LENGTH - 1)
    {
      forward_name[NAME_STR_LENGTH - 1] = '\0';
    }
    status = sprintf(reverse_name, "%d_r", reaction_definition_count);
    if (status < 0) return FALSE;
    if (status > NAME_STR_LENGTH - 1)
    {
      reverse_name[NAME_STR_LENGTH - 1] = '\0';
    }
#endif
  }
  else
  {
    name_len = strlen(name) + 1;
    // Function sizeof() includes terminating NUL in count.
    forward_name = malloc(name_len + sizeof(FORWARD_SUFFIX));
    strcpy(forward_name, name);
    strcat(forward_name, FORWARD_SUFFIX);
    reverse_name = malloc(name_len + sizeof(REVERSE_SUFFIX));
    strcpy(reverse_name, name);
    strcat(reverse_name, REVERSE_SUFFIX);
  }

  if (1 == extra_rhs_molecules->count)
  {
    rule_t *new_molecule_rule;
    int new_specie_id;

    // Need to change to number of species later.
    rhs_count -= extra_rhs_molecules->count;

    new_molecule_rule =
      world->dynrules->rules.a + extra_rhs_molecules->ints.a[0];

    status =
      rule_create_specie(new_molecule_rule, &new_specie_id, world->dynspecies,
                         world->dynrules);
    if (0 == status)
    {
      return FALSE;
    }
    DEBUG_PRINTF(2, "Rule %s creates new specie: ", forward_name);
    CGRAPH_OUTPUT_BNGL(OUTPUT_DEBUG, 2,
                       (world->dynspecies->species.a + new_specie_id)->cg,
                       world->dyncolors);
    DEBUG_PRINTF(2, "\n");

    /* Keep track of the species that need to be created. */
    status = dynints_int_add(new_specie_id, extra_rhs_species);
    if (0 == status)
    {
      return FALSE;
    }
  }


  dynints_destroy(extra_lhs_molecules);
  dynints_destroy(extra_rhs_molecules);

  if ((0 == extra_lhs_species->count) ||
      (FALSE == bidirectional))
  {
    dynints_destroy(extra_lhs_species);
    extra_lhs_species = NULL;
  }
  if (0 == extra_rhs_species->count)
  {
    dynints_destroy(extra_rhs_species);
    extra_rhs_species = NULL;
  }

#ifdef RULE_CREATE_DEBUG
  DYNRULES_OUTPUT(OUTPUT_DEBUG, 6, world->dynrules, world->dyncolors);
#endif

  /* Setup forward reaction. */
  status = setup_reaction_rule(forward_name, forward_rate, flags, lhs_count,
                               rhs_count, specie_rule_indexes,
                               extra_rhs_species);
  if (FALSE == status)
  {
    yyerror_prev("Error setting up forward reaction");
    dynints_destroy(extra_lhs_species);
    dynints_destroy(extra_rhs_species);
    return FALSE;
  }

  if (TRUE == bidirectional)
  {
    /* Setup reverse reaction. */
    status = setup_reaction_rule(reverse_name, reverse_rate, flags, rhs_count,
                                 lhs_count, specie_rule_indexes + lhs_count,
                                 extra_lhs_species);
    if (FALSE == status)
    {
      yyerror_prev("Error setting up reverse reaction");
    dynints_destroy(extra_rhs_species);
      return FALSE;
    }
  }

  if (NULL != forward_name)
  {
    free(forward_name);
  }
  if (NULL != reverse_name)
  {
    free(reverse_name);
  }
  return TRUE;
}

boolean
finish_observable(const char *name, const int rule_index, double denominator)
{
  int status;

#ifdef RULE_CREATE_DEBUG
  DYNRULES_OUTPUT(OUTPUT_DEBUG, 10, world->dynrules, world->dyncolors);
  DYNCONNECTIONS_OUTPUT(OUTPUT_DEBUG, 3, dynconnections);
#endif

  // Applies connections between molecules.
  status = dynconnections_apply_to_rules(1, 0, ds_system.strict,
                                         specie_rule_indexes, world->dyncolors,
                                         world->dynrules, dynconnections);
  if (FALSE == status)
  {
    yyerror_prev("Error applying connections to observable rules");
    return FALSE;
  }

  status = report_definition_add(name, rule_index, denominator, world);

  return status;
}

boolean
setup_reaction_rule(const char *name, const double rate, const int flags,
                    const int reactant_count, const int product_count,
                    int *current_specie_rule_indexes,
                    dynints_t *new_species)
{
  int status;
  reactant_definition_t *reactant_definitions;
  int i;

  status =
    reaction_definition_add(name, rate, flags, reactant_count,
                            product_count, new_species,
                            &world->reaction_definition_list);
  if (FALSE == status)
  {
    yyerror_prev("Error normalizing specie");
    return FALSE;
  }

  /* New reaction is at head of linked list. */
  reactant_definitions = world->reaction_definition_list->reactant_definitions;

  /* Add reactant rules. */
  for (i = 0; i < reactant_count; i++)
  {
    reactant_definitions[i].rule_index = current_specie_rule_indexes[i];
  }

  return TRUE;
}

boolean
bng_command_set_default()
{
  if (NULL != bng_command_args.prefix)
  {
    free(bng_command_args.prefix);
    bng_command_args.prefix = NULL;
  }
  if (NULL != bng_command_args.suffix)
  {
    free(bng_command_args.suffix);
    bng_command_args.suffix = NULL;
  }
  bng_command_args.prefix = NULL;
  bng_command_args.atol = 1e-8;
  bng_command_args.rtol = 1e-8;
  bng_command_args.sparse = 0;
  bng_command_args.n_steps = 1.0;
  bng_command_args.t_start = 0.0;
  bng_command_args.t_end = 0.0;
  bng_command_args.seed = 0;
  bng_command_args.volume = 0.0;
  bng_command_args.time_delta = 0;
  bng_command_args.max_prob = 0.999999;
  bng_command_args.pseudo_factor = 0.0;
  bng_command_args.pseudo_count = 0.0;
  bng_command_args.suffix = NULL;

  return TRUE;
}

boolean
bng_command_set_double(const char *name, const double value)
{
  if (0 == strcasecmp("atol", name))
  {
    bng_command_args.atol = value;
  }
  else if (0 == strcasecmp("rtol", name))
  {
    bng_command_args.rtol = value;
  }
  else if (0 == strcasecmp("sparse", name))
  {
    bng_command_args.sparse = (int)value;
  }
  else if (0 == strcasecmp("n_steps", name))
  {
    bng_command_args.n_steps = value;
  }
  else if (0 == strcasecmp("t_start", name))
  {
    bng_command_args.t_start = value;
  }
  else if (0 == strcasecmp("t_end", name))
  {
    bng_command_args.t_end = value;
  }
  else if (0 == strcasecmp("seed", name))
  {
    bng_command_args.seed = (unsigned)value;
  }
  else if (0 == strcasecmp("volume", name))
  {
    bng_command_args.volume = value;
  }
  else if (0 == strcasecmp("time_delta", name))
  {
    bng_command_args.time_delta = value;
  }
  else if (0 == strcasecmp("max_prob", name))
  {
    bng_command_args.max_prob = value;
  }
  else if (0 == strcasecmp("pseudo_factor", name))
  {
    bng_command_args.pseudo_factor = value;
  }
  else if (0 == strcasecmp("pseudo_count", name))
  {
    bng_command_args.pseudo_count = value;
  }
  else
  {
    error_printf("Unknown floating point parameter: %s\n", name);
    yyerror("Unknown parameter");
    return FALSE;
  }

  return TRUE;
}

boolean
bng_command_set_string(const char *name, const char *value)
{
  if (0 == strcasecmp("prefix", name))
  {
    bng_command_args.prefix = strdup(value);
    error_printf("Unknown string parameter: %s\n", name);
    yyerror("Unknown parameter");
    return FALSE;
  }
  else if (0 == strcasecmp("suffix", name))
  {
    bng_command_args.suffix = strdup(value);
    error_printf("Unknown string parameter: %s\n", name);
    yyerror("Unknown parameter");
    return FALSE;
  }
  else
  {
    error_printf("Unknown string parameter: %s\n", name);
    yyerror("Unknown parameter");
    return FALSE;
  }

  return TRUE;
}

boolean
setup_ps(bng_command_args_t *args)
{
  int status;

  status = world_complete(world);

  return status;
}

boolean
simulate_rs(bng_command_args_t *args)
{
  int status;
  double logging_interval;

  MYSRAND(args->seed);

  if (args->n_steps <= 0)
  {
    // Prints a log for every generation.
    logging_interval = 0.0;
  }
  else
  {
    //RG: Use world->current_time rather than args->t_start to set
    //    loggging_interval. It appears t_start is not used in any other way,
    //    so I think this change is safe.
    logging_interval = (args->t_end - world->current_time) / (args->n_steps);
  }

  if(args->t_start > 0)
  {
    //RG: Warn that t_start is now ignored.
    printf("Warning, t_start option is ignored by simulate_rs.\n");
  }

  debug_printf(2,
             "In simulate_rs, t_end:%e, volume:%e, time_delta:%e, "
             "max_prob:%e\n",
             args->t_end, args->volume, args->time_delta, args->max_prob);

  if (NULL != world->report_definition_list)
  {
    status = reactionstoc_run(world, logging_interval,
                              report_definition_output_all_reports,
                              args->t_end,
                              args->pseudo_factor, args->pseudo_count);
  }
  else
  {
    status = reactionstoc_run(world, logging_interval,
                              world_output_pop, args->t_end,
                              args->pseudo_factor, args->pseudo_count);
  }

  return status;
}

boolean
simulate_ps(bng_command_args_t *args)
{
  int status;
  double logging_interval;

  MYSRAND(args->seed);

  if (args->n_steps <= 0)
  {
    // Prints a log for every generation.
    logging_interval = 0.0;
  }
  else
  {
    logging_interval = (args->t_end - args->t_start) / (args->n_steps);
  }

  debug_printf(2,
             "In simulate_ps, t_end:%e, volume:%e, time_delta:%e, "
             "max_prob:%e\n",
             args->t_end, args->volume, args->time_delta, args->max_prob);

  if (NULL != world->report_definition_list)
  {
    status = particlestoc_run(world, logging_interval,
                              report_definition_output_all_reports,
                              args->t_end, args->volume,
                              args->time_delta, args->max_prob,
                              args->pseudo_factor,
                              args->pseudo_count, ds_system.avogadro);
  }
  else
  {
    status = particlestoc_run(world, logging_interval,
                              world_output_pop,
                              args->t_end, args->volume,
                              args->time_delta, args->max_prob,
                              args->pseudo_factor,
                              args->pseudo_count, ds_system.avogadro);
  }


  return status;
}

boolean
simulate_ss(bng_command_args_t *args)
{
  int status;
  double logging_interval;

  status = FALSE;

  MYSRAND(args->seed);

  if (0 == args->n_steps)
  {
    // Prints a log for every generation.
    logging_interval = 0.0;
  }
  else
  {
    logging_interval = (args->t_end - args->t_start) / args->n_steps;
  }

  if (NULL != world->report_definition_list)
  {
#ifdef ENABLE_SPECIESTOC
    status = dynstoc(world, logging_interval,
                     report_definition_output_all_reports,
                     args->t_end);
#endif /* ENABLE_SPECIESTOC */
  }
  else
  {
#ifdef ENABLE_SPECIESTOC
    status = dynstoc(world, logging_interval,
                     world_output_pop, args->t_end);
#endif /* ENABLE_SPECIESTOC */
  }


  return status;
}

boolean
state_output_ps(bng_command_args_t *args)
{
  particlestoc_output_all(OUTPUT_STANDARD, 0, world, world->current_time);

  return TRUE;
}

static void ds_system_setup()
{
  // Default to strict syntax checking.
  ds_system.strict = 1;

  ds_system.avogadro = AVOGADRO;
}

boolean
ds_system_set_value(const char *name, double value)
{
  if (0 == strcasecmp(name, "verbose"))
  {
    set_debug_level(value);
  }
  else if (0 == strcasecmp(name, "strict"))
  {
    ds_system.strict = value;
  }
  else if (0 == strcasecmp(name, "avogadro"))
  {
    ds_system.avogadro = value;
  }
  else
  {
    return FALSE;
  }

  return TRUE;
}

boolean
ds_system_get_value(const char *name, double *pvalue)
{
  if (0 == strcasecmp(name, "avogadro"))
  {
    *pvalue = ds_system.avogadro;
  }
  else
  {
    return FALSE;
  }

  return TRUE;
}
