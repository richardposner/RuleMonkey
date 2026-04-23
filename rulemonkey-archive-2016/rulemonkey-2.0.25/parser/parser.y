
%{
#include <math.h>
#include <string.h>
#include <stdio.h> /* For YYDEBUG. */
#include <assert.h>

#include "parser_lib.h"
#include "constants.h"
#include "dynsymbols.h"
#include "dyncolors.h"
#include "dynmolecules.h"
#include "dynconnections.h"
#include "dynmaps.h"
#include "dynspecies.h"
#include "cgraph.h"
#include "reactant_definition.h"
#include "world.h"
#include "report_definition.h"
#include "output.h"

#define YYERROR_VERBOSE

int yylex(void); // To avoid gcc warnings.

void suppress_warnings29() {SUPPRESS_WARNINGS;}

#if ac_YYDEBUG
#define YYDEBUG 1
#endif

%}

%union {
  int ival;
  double dval;
  char str[1024];
}
%token <str> LABEL
%token <str> QSTRING
%token <str> NAME
%token <str> NUMBER
%token QUESTION_MARK
%token MOLECULES SPECIES
%token BEGIN_PARAMETERS END_PARAMETERS
%token BEGIN_SPECIES END_SPECIES
%token BEGIN_REACTION_RULES END_REACTION_RULES
%token BEGIN_OBSERVABLES END_OBSERVABLES
%token BEGIN_MOLECULE_TYPES END_MOLECULE_TYPES
%token BEGIN_ACTIONS END_ACTIONS
%token UNIDIRECTIONAL_RULE BIDIRECTIONAL_RULE
%token GENERATE_NETWORK SIMULATE_ODE SIMULATE_SSA TO_SBML_FILE
%token SETUP_PS SIMULATE_PS SIMULATE_RS STATE_OUTPUT_PS
%token CMP_LE CMP_GE CMP_LT CMP_GT
%token LOGICAL_AND LOGICAL_OR LOGICAL_XOR
%left '-' '+'
%left '*' '/'
%left UMINUS  /* negation--unary minus */
%right '^' /* exponentiation */

%type <ival> molecule_component_name
%type <ival> molecule_component_modifier_list molecule_component_modifier
%type <ival> optional_molecule_component_modifier_list
%type <ival> optional_object_modifier_list
%type <ival> object_modifier_list object_modifier molecule_link_modifier
%type <ival> specie_rule_plus_list
%type <ival> reaction_rule_lhs reaction_rule_rhs

%type <dval> expression expression_element

%type <str> symbol
%type <str> name_or_number
%type <str> name_or_number_or_wildcard
%type <str> optional_label_and_or_number
%type <ival> molecule_name
%type <ival> molecule_name_declaration
%type <ival> molecule_modified
%type <ival> molecule
%type <ival> molecule_declaration
%type <ival> molecule_aggregate
%type <ival> molecule_aggregate_declaration
%type <ival> molecule_aggregate_modified
%type <ival> molecule_aggregate_modified_declaration
%type <ival> specie_declaration
%type <ival> specie
%type <ival> specie_rule
%type <ival> specie_rule_or_list
%type <str> bng_command_name

%%

bngl_file:
    bng_item_list
      { ; }
  ;

bng_item_list:
    bng_item
      { ; }
  | bng_item_list bng_item
      { ; }
  ;

bng_item:
    bng_command
      { ; }
  | parameter_block
      { ; }
  | specie_block
      { ; }
  | reaction_rule_block
      { ; }
  | observable_block
      { ; }
  | molecule_types_block
      { ; }
  | action_block
      { ; }
  ;

bng_command:
    bng_command_name '(' optional_perl_hash ')' ';'
      {
        boolean status;

        if (*$1 != '\0')
        {
          /* Execute command, passing collected parameters. */
          status = bng_command_ptr(&bng_command_args);
          if (FALSE == status)
          {
            yyerror("Error running bng command");
            exit(1);
          }
        }
      }
  ;

optional_perl_hash:
    /* Empty. */
  | perl_hash
      { ; }
  ;

bng_command_name:
    GENERATE_NETWORK
      {
        yyerror_prev("generate_network() not implemented.");
        *$$ = '\0';
      }
  | SIMULATE_ODE
      {
        yyerror_prev("simulate_ode() not implemented.");
        *$$ = '\0';
      }
  | SIMULATE_SSA
      {
        yyerror_prev("simulate_ssa() not implemented.");
        *$$ = '\0';
      }
  | STATE_OUTPUT_PS
      {
        /* Set callback. */
        bng_command_ptr = state_output_ps;
        bng_command_set_default();
        strcpy($$, "ps_state_output");
      }
  | SETUP_PS
      {
        /* Set callback. */
        bng_command_ptr = setup_ps;
        bng_command_set_default();
        strcpy($$, "setup_ps");
      }
  | SIMULATE_PS
      {
// simulate_ps does not currently work, rule checking has infinite loop.
#ifdef ENABLE_SIMULATE_PS
        // simulate_ps does not currently work, rule checking has infinite loop.
        /* Set callback. */
        bng_command_ptr = simulate_ps;
        bng_command_set_default();
        strcpy($$, "simulate_ps");
#else
        // simulate_ps does not currently work, rule checking has infinite loop.
        yyerror_prev("simulate_ps() not implemented.");
        *$$ = '\0';
#endif
      }
  | SIMULATE_RS
      {
        /* Set callback. */
        bng_command_ptr = simulate_rs;
        bng_command_set_default();
        strcpy($$, "simulate_rs");
      }
  | TO_SBML_FILE
      {
        yyerror_prev("to_sbml_file() not implemented.");
        *$$ = '\0';
      }
  ;

parameter_block:
    begin_parameters end_parameters
      { ; }
  | begin_parameters parameter_block_list end_parameters
      { ; }
  ;

specie_block:
    begin_species end_species
      { ; }
  | begin_species specie_block_list end_species
      { ; }
  ;

reaction_rule_block:
    begin_reaction_rules end_reaction_rules
      { ; }
  | begin_reaction_rules reaction_rule_block_list end_reaction_rules
      { ; }
  ;

observable_block:
    begin_observables end_observables
      { ; }
  | begin_observables observable_block_list end_observables
      { ; }
  ;

molecule_types_block:
    begin_molecule_types end_molecule_types
      { ; }
  | begin_molecule_types molecule_type_block_list end_molecule_types
      { ; }
  ;

action_block:
    begin_actions end_actions
      { ; }
  | begin_actions action_block_list end_actions
      { ; }
  ;

begin_parameters:
    BEGIN_PARAMETERS
      { current_block = PARAMETERS_BLOCK_TYPE; }
    ;

end_parameters:
    END_PARAMETERS
      { current_block = NO_BLOCK_TYPE; }
    ;

begin_species:
    BEGIN_SPECIES
      { current_block = SPECIES_BLOCK_TYPE; }
    ;

end_species:
    END_SPECIES
      { current_block = NO_BLOCK_TYPE; }
    ;

begin_reaction_rules:
    BEGIN_REACTION_RULES
      { current_block = REACTION_RULES_BLOCK_TYPE; }
    ;

end_reaction_rules:
    END_REACTION_RULES
      { current_block = NO_BLOCK_TYPE; }
    ;

begin_observables:
    BEGIN_OBSERVABLES
      { current_block = OBSERVABLES_BLOCK_TYPE; }
    ;

end_observables:
    END_OBSERVABLES
      { current_block = NO_BLOCK_TYPE; }
    ;

begin_molecule_types:
    BEGIN_MOLECULE_TYPES
      {
        current_block = MOLECULE_TYPES_BLOCK_TYPE;
      }
    ;

end_molecule_types:
    END_MOLECULE_TYPES
      {
        current_block = NO_BLOCK_TYPE;
        molecule_types_present = 1;
      }
    ;

begin_actions:
    BEGIN_ACTIONS
      { current_block = ACTIONS_BLOCK_TYPE; }
    ;

end_actions:
    END_ACTIONS
      { current_block = NO_BLOCK_TYPE; }
    ;

parameter_block_list:
    parameter_block_list_item
      { ; }
  | parameter_block_list parameter_block_list_item
      { ; }
  ;

specie_block_list:
    specie_block_list_item
      { ; }
  | specie_block_list specie_block_list_item
      { ; }
  ;

reaction_rule_block_list:
    reaction_rule_block_list_item
      { ; }
  | reaction_rule_block_list reaction_rule_block_list_item
      { ; }
  ;

observable_block_list:
    observable_block_list_item
      { ; }
  | observable_block_list observable_block_list_item
      { ; }
  ;

molecule_type_block_list:
    molecule_type_block_list_item
      { ; }
  | molecule_type_block_list molecule_type_block_list_item
      { ; }
  ;

action_block_list:
    action_block_list_item
      { ; }
  | action_block_list action_block_list_item
      { ; }
  ;

parameter_block_list_item:
    optional_label_and_or_number parameter_declaration
      {
      }
  ;

parameter_declaration:
    symbol expression
      {
        int status;

        if (0 == strncasecmp($1, SYSTEM_PREFIX, sizeof(SYSTEM_PREFIX) - 1))
        {
          char *system_symbol;

          system_symbol = $1;

          // Skip over prefix.
          system_symbol += sizeof(SYSTEM_PREFIX) - 1;

          status = ds_system_set_value(system_symbol, $2);
          if (0 == status)
          {
            yyerror("System value specified on previous line is not valid");
            exit(1);
          }
        }
        else
        {
          status = dynsymbols_name_update($1, $2, dynsymbols);
          if (0 == status)
          {
            yyerror("Error creating new parameter");
            assert(0);
            exit(1);
          }
        }
      }
  | symbol '=' expression
      {
        int status;

        if (0 == strncasecmp($1, SYSTEM_PREFIX, sizeof(SYSTEM_PREFIX) - 1))
        {
          char *system_symbol;

          system_symbol = $1;

          // Skip over prefix.
          system_symbol += sizeof(SYSTEM_PREFIX) - 1;

          status = ds_system_set_value(system_symbol, $3);
          if (0 == status)
          {
            yyerror("System value specified on previous line is not valid");
            exit(1);
          }
        }
        else
        {
          status = dynsymbols_name_update($1, $3, dynsymbols);
          if (0 == status)
          {
            yyerror("Error creating new parameter");
            assert(0);
            exit(1);
          }
        }
      }
  ;

symbol:
    NAME
  ;

action_block_list_item:
    bng_command
      { ; }
  ;

specie_block_list_item:
    optional_label_and_or_number specie_declaration expression
      {
        boolean status;
        boolean new;
        int specie_id;
        double rounded_pop;

        /* Finished initializing specie. */
        initializing_specie = FALSE;

        /* Add inter-molecule connections to cgraph. */
        dynconnections_apply_to_graph
          (current_specie_cgraph, dynconnections);

#ifdef HAVE_ROUND
        rounded_pop = round($3);
#else
        rounded_pop = floor($3 + 0.5);
#endif

#ifndef CREATE_NON_AGGREGATE_SPECIE
        status = dynspecies_specie_add(current_specie_cgraph, $1, rounded_pop,
                                      TRUE, TRUE, FALSE, FALSE,
                                      &specie_id, &new, world->dynspecies,
                                      world);
        if (FALSE == status)
        {
          yyerror_prev("Error adding specie");
          assert(0);
          exit(1);
        }
#else
        {
          size_t i;
          size_t specie_index;

          specie_index = rounded_pop;
          if (specie_index != rounded_pop)
          {
            yyerror("New population too large");
            assert(0);
            exit(1);
          }


          for (i = 0; i < specie_index; i++)
          {
            status = dynspecies_specie_add(current_specie_cgraph, $1,
                                           rounded_pop,
                                           TRUE, TRUE, FALSE, FALSE,
                                           &specie_id, &new, world->dynspecies,
                                           world);
            if (FALSE == status)
            {
              yyerror_prev("Error adding specie");
              assert(0);
              exit(1);
            }
          }
        }
#endif
      }
  ;

molecule_type_block_list_item:
    molecule_declaration
      {
        /* Finished initializing specie (always a single molecule). */
        initializing_specie = FALSE;

        /* Single molecule, so no inter-molecule connection to add. */

        /* Just a molecule declaration, so do not add to specie list. */
      }
  ;
     

reaction_rule_block_list_item:
    reaction_rule
      { ; }
  ;

reaction_rule:
    reaction_rule_type
      {
        /* Finished building rule. */
        building_rule = FALSE;
      }
  ;

reaction_rule_type:
    optional_label_and_or_number reaction_rule_lhs UNIDIRECTIONAL_RULE reaction_rule_rhs expression
      {
        int status;

        status = finish_reaction_rules($1, FALSE, $2, $4, $5, 0.0);
        if (FALSE == status)
        {
          yyerror_prev("Error setting up reaction rule");
          assert(0);
          exit(1);
        }

      }
  | optional_label_and_or_number reaction_rule_lhs BIDIRECTIONAL_RULE reaction_rule_rhs expression ',' expression
      {
        int status;

        status = finish_reaction_rules($1, TRUE, $2, $4, $5, $7);
        if (FALSE == status)
        {
          yyerror_prev("Error setting up reaction rule");
          exit(1);
        }
      }
  ;

reaction_rule_lhs:
    specie_rule_plus_list
      {
        /* Finished with left hand side. */
        reaction_lhs = FALSE;

        molecule_star_lhs = molecule_star_count;
        molecule_star_count = 0;

        /* Return size of list. */
        $$ = $1;
      }
  ;

reaction_rule_rhs:
    specie_rule_plus_list
      {
        /* Finished with right hand side. */

        molecule_star_rhs = molecule_star_count;
        molecule_star_count = 0;

        /* Return size of list. */
        $$ = $1;
      }
  ;

specie_rule_or_list:
    specie_rule
      {
        /* Keep track of connections for each specie. */
        building_rule = TRUE;

        $$ = $1;
      }
  | specie_rule_or_list specie_rule
      {
        /* Make rule 'or' to lhs rule. */
        dynrules_rule_add_or_rule($1, $2, world->dynrules);

        $$ = $1;
      }
  ;

specie_rule_plus_list:
    '*'
      {
        /* Specie is getting created or deleted. */
        molecule_star_count++;

        $$ = 0;
      }
  | specie_rule
      {
        /* Keep track of connections for each specie. */
        building_rule = TRUE;

        /* Return size of list. */
        $$ = 1;
      }
  | specie_rule_plus_list '+' '*'
      {
        /* Specie is getting created or deleted. */
        molecule_star_count++;

        /* Return size of list. */
        $$ = $1;
      }
  | specie_rule_plus_list '+' specie_rule
      {
        /* Return size of list. */
        $$ = $1 + 1;
      }
  ;

observable_block_list_item:
    optional_label_and_or_number MOLECULES NAME specie_rule_or_list
      {
        int status;

        building_rule = FALSE;

        // Passing a denominator of one means count each molecule.
        status = finish_observable($3, $4, 1.0);
        if (FALSE == status)
        {
          yyerror_prev("Finishing molecule observable failed.");
          assert(0);
          exit(1);
        }
      }
  | optional_label_and_or_number SPECIES NAME specie_rule_or_list
      {
        int status;

        building_rule = FALSE;

        // Passing a denominator of zero means just count species,
        // no matter how many matching molecules found in each species.
        status = finish_observable($3, $4, 0.0);
        if (FALSE == status)
        {
          yyerror_prev("Finishing specie observable failed.");
          assert(0);
          exit(1);
        }
      }
  ;

specie_rule:
    '(' specie_rule ')'
      {
        // Return last molecule rule index.
        $$ = $2;
      }
  | specie
      {
        // Return last molecule rule index.
        $$ = $1;
      }
  | specie_rule LOGICAL_AND specie
      {
        /* Make rule 'and' to lhs rule. */
        dynrules_rule_add_and_rule($1, $3, world->dynrules);

        $$ = $1;
      }
  | specie_rule LOGICAL_OR specie
      {
        /* Make rule 'or' to lhs rule. */
        dynrules_rule_add_or_rule($1, $3, world->dynrules);

        // Return last molecule rule index.
        $$ = $1;
      }
  | specie_rule LOGICAL_XOR specie
      {
        /* Make rule 'xor' to lhs rule. */
        dynrules_rule_add_xor_rule($1, $3, world->dynrules);

        // Return last molecule rule index.
        $$ = $1;
      }
  ;

specie:
    molecule_aggregate_modified
      {
        /* Finished building current specie. */
        building_specie = FALSE;

        // Return last molecule index.
        $$ = $1;
      }
  ;

specie_declaration:
    molecule_aggregate_modified_declaration
      {
        /* Finished building current specie. */
        building_specie = FALSE;

        // Return last molecule index.
        $$ = $1;
      }
  ;

molecule_aggregate_modified:
    molecule_aggregate
      {
        $$ = $1;
      }
/*
  | '(' molecule_aggregate ')' object_modifier_list
      {
        yyerror_prev("Modified molecule aggregate not yet supported");
        assert(0);
        exit(1);

        $$ = $2;
      }
*/
  ;

molecule_aggregate_modified_declaration:
    molecule_aggregate_declaration
      {
        $$ = $1;
      }
  | '(' molecule_aggregate_declaration ')' optional_object_modifier_list
      {
        yyerror_prev("Modified molecule aggregate not yet supported");
        assert(0);
        exit(1);

        $$ = $2;
      }
  ;

molecule_aggregate_declaration:
    molecule_declaration
      {
        /* Return lhs molecule index. */
        $$ = $1;
      }
  | molecule_aggregate_declaration '.' molecule_declaration
      {
        // Don't add_and_rule because it breaks specie declaration for species with multiple molecules.
        //dynrules_rule_add_and_rule($1, $3, world->dynrules);
        /* Return lhs molecule index. */
        $$ = $1;
      }
  ;

molecule_aggregate:
    molecule_modified
      {
        /* Return lhs molecule index. */
        $$ = $1;
      }
  | molecule_aggregate '.' molecule_modified
      {
        dynrules_rule_add_and_rule($1, $3, world->dynrules);

        /* Return lhs molecule index. */
        $$ = $1;
      }
  ;

molecule_modified:
    molecule optional_object_modifier_list
      {
        /* Finished with current molecule. */
        current_molecule = NULL;
        current_molecule_index = -1;
        current_molecule_rule_index = -1;
        current_object_index = -1;
        current_object_rule_index = -1;
        new_molecule = FALSE;

        /* Return molecule index. */
        $$ = $1;
      }
  ;

molecule:
    molecule_name molecule_components
      {
        free(molecule_component_multipliers_table);
        molecule_component_multipliers_table = NULL;

        /* Return molecule index. */
        $$ = $1;
      }
  ;

molecule_declaration:
    molecule_name_declaration molecule_components
      {
        free(molecule_component_multipliers_table);
        molecule_component_multipliers_table = NULL;

        /* Return molecule index. */
        $$ = $1;
      }
  ;

molecule_name:
    NAME
      {
        int molecule_color;
        boolean status;

        /* Create molecule to verify structure. */
        molecule_color = setup_molecule_name($1, FALSE);

        if (FALSE == building_specie)
        {
          if (FALSE == building_rule)
          {
            /* This might be the first molecule in a reaction rule. */

            /* Reset inter-molecule connections. */
            dynconnections_reset(dynconnections);

            /* Reset reaction-product mappings. */
            dynmaps->count = 0;

            /* Reset specie counters. */
            current_specie_index = 0;
            reaction_lhs = TRUE;
          }
          else
          {
            // ((TRUE == building_rule) && (FALSE == building_specie))
            // Starting to build the 2nd or 3rd specie (in a reaction rule)
            current_specie_index++;
          }
        }

        /* For creating rules. */
        current_molecule_name = $1;

        /* Add rule to match molecule with the given name. */
        /* Count may be modified later. */
        /* count >=1 to allow rules to be multiplexed by default. */
        status = dynrules_rule_create(molecule_color, NO_MODIFIER_FLAGS,
                                      GREAT_EQ, 1, 1.0, -1, -1, -1, -1,
                                      world->dynrules);
        if (FALSE == status)
        {
            yyerror_prev("Unable to create rule for molecule");
            assert(0);
            exit(1);
        }

        /* Save molecule rule index. */
        current_molecule_rule_index = world->dynrules->count - 1;

        if (FALSE == building_specie)
        {
          /* Beginning of new specie. */
          building_specie = TRUE;

          /* Save master rule just created for new specie. */
          if (current_specie_index >= MAX_SPECIE_RULE_INDEX)
          {
            yyerror_prev("Too many species on same rule, increase "
                         "MAX_SPECIE_RULE_INDEX");
            assert(0);
            exit(1);
          }

          specie_rule_indexes[current_specie_index] =
            current_molecule_rule_index;
        }

        /* Return molecule index for use in building rules. */
        $$ = current_molecule_rule_index;
      }
  ;

molecule_name_declaration:
    NAME
      {
        int color;

        /* Create molecule to define structure or */
        /* verify structure if molecule already exists. */
        color = setup_molecule_name($1, TRUE);

        if (FALSE == initializing_specie)
        {
          /* Beginning of new specie initialization. */
          initializing_specie = TRUE;

          /* Initialize molecule size to zero. */
          cgraph_reset(current_specie_cgraph);

          /* Reset inter-molecule connections. */
          dynconnections_reset(dynconnections);
          current_specie_index = 0;
          reaction_lhs = 0;
        }

        /* Add primary molecule vertex to graph. */
        current_molecule_index =
          cgraph_add_node(current_specie_cgraph, color);

        /* Return molecule index for use in building rules. */
        $$ = current_molecule_index;
      }
  ;

molecule_components:
    '(' molecule_component_list ')'
      {
        /*
         * Now that components are finished, molecule modifiers will be
         * read if present.
         */
        current_object_index = current_molecule_index;
        current_object_name = current_molecule_name;
        current_object_rule_index = current_molecule_rule_index;
      }
  ;

molecule_component_list:
    /* Empty */
  | molecule_component
      { ; }
  | molecule_component_list ',' molecule_component
      { ; }
  ;

molecule_component:
    molecule_component_name optional_molecule_component_modifier_list
      {
        boolean status;
        dynints_t *state_dynints;

        if (TRUE == initializing_specie)
        {
          state_dynints = current_specie_cgraph->
                          state_dynints_array.a[current_object_index];

          if (TRUE == new_molecule)
          {
            status = molecule_component_add($1, $2, state_dynints,
                                            world->dyncolors, current_molecule);
            if (FALSE == status)
            {
              yyerror("Unable to add molecule_component");
              exit(1);
            }
          }
          else
          {
            status =
              molecule_component_verify($1, $2, state_dynints,
                                        molecule_component_multipliers_table,
                                        world->dyncolors, current_molecule,
                                        molecule_types_present);
            if (FALSE == status)
            {
              yyerror("Unable to verify molecule_component");
              exit(1);
            }
          }
        }
        else
        {
          rule_t *object_rule;

          object_rule = world->dynrules->rules.a + current_object_rule_index;

          status =
            molecule_component_verify($1, $2,
                                      object_rule->state_dynints,
                                      molecule_component_multipliers_table,
                                      world->dyncolors, current_molecule,
                                      molecule_types_present);
          if (FALSE == status)
          {
            yyerror("Unable to verify molecule_component");
            exit(1);
          }
        }
      }
  ;

molecule_component_name:
    NAME
      {
        color_t *component_color;

        /* Lookup component color. */
        component_color = real_dyncolors_color_create($1, COLOR_TYPE_COMPONENT,
                                                      world->dynspecies,
                                                      world->dyncolors, world);
        if (NULL == component_color)
        {
          yyerror_prev("Error initializing component color");
          assert(0);
          exit(1);
        }

        if (TRUE == initializing_specie)
        {
          /* Add component node. */
          current_component_index =
            cgraph_add_node(current_specie_cgraph, component_color->value);
          if (-1 == current_component_index)
          {
            yyerror_prev("Unable to add component to specie");
            assert(0);
            exit(1);
          }

          /* For assigning object modifiers. */
          current_object_index = current_component_index;

          /* Connect component to molecule. */
          cgraph_connect_nodes(current_specie_cgraph, current_component_index,
                               current_molecule_index);
        }
        else
        {
          int status;

          /* Building rules to match specie specification. */
          /* Count may be modified later. */
          /* count >=1 to allow rules to be multiplexed by default. */
          status = dynrules_rule_create(component_color->value,
                                        ADJACENT_EXCLUSIVE_MODIFIER,
                                        GREAT_EQ, 1, 1.0, -1, -1, -1, -1,
                                        world->dynrules);
          if (FALSE == status)
          {
              yyerror_prev("Unable to create rule for molecule");
              assert(0);
              exit(1);
          }

          /* So that object modifiers know what they are modifying. */
          current_object_rule_index = world->dynrules->count - 1;
          current_object_name = $1;

          /* Associate component with molecule. */
          dynrules_rule_add_adjacent_rule (current_molecule_rule_index,
                                           current_object_rule_index,
                                           world->dynrules);
        }

        $$ = component_color->value;
      }
  ;

optional_molecule_component_modifier_list:
    /* Empty */
      {
        $$ = 1;
      }
  | molecule_component_modifier_list
      {
        $$ = $1;
      }
  ;

molecule_component_modifier_list:
    molecule_component_modifier
      {
        if (-1 == $1)
        {
          $$ = 1;
        }
        else
        {
          $$ = $1;
        }
      }
  | molecule_component_modifier_list molecule_component_modifier
      {
        if (-1 == $1)
        {
          $$ = $2;
        }
        else if (-1 == $2)
        {
          $$ = $1;
        }
        else
        {
          /* Two valid component multiples, so multiply them together. */
          $$ = $1 * $2;
        }
      }
  ;

optional_object_modifier_list:
    /* Empty */
      {
        $$ = 1;
      }
  | object_modifier_list
      {
        $$ = $1;
      }
  ;

object_modifier_list:
    object_modifier
      {
        if (-1 == $1)
        {
          $$ = 1;
        }
        else
        {
          $$ = $1;
        }
      }
  | object_modifier_list object_modifier
      {
        if (-1 == $1)
        {
          $$ = $2;
        }
        else if (-1 == $2)
        {
          $$ = $1;
        }
        else
        {
          /* Two valid component multiples, so multiply them together. */
          $$ = $1 * $2;
        }
      }
  ;

molecule_link_modifier:
    '+'
      {
        /* Add link wildcard. */
        if (initializing_specie)
        {
          yyerror_prev("Link wildcard not allowed during specie initialization");
          assert(0);
          exit(1);
        }
        else
        {
          /* Must match one or more bonds. */
          dynrules_rule_create(COLOR_WILDCARD, NO_MODIFIER_FLAGS, GREAT_EQ, 1,
                               1.0, -1, -1, -1, -1, world->dynrules);
          dynrules_rule_add_adjacent_rule
                    (current_object_rule_index,
                     (int)world->dynrules->count - 1,
                     world->dynrules);
        }

        $$ = -1;
      }
  | QUESTION_MARK
      {
        /* Add link wildcard. */
        if (initializing_specie)
        {
          yyerror_prev("Link wildcard not allowed during specie initialization");
          assert(0);
          exit(1);
        }
        else
        {
          /* Must match zero or one bond. */
          dynrules_rule_create(COLOR_WILDCARD, NO_MODIFIER_FLAGS, LESS_EQ, 1,
                               1.0, -1, -1, -1, -1, world->dynrules);
          dynrules_rule_add_adjacent_rule
                     (current_object_rule_index,
                      (int)world->dynrules->count - 1,
                      world->dynrules);
        }

        $$ = -1;
      }
  | name_or_number
      {
        int status;

        /* Add link name to component. */
        if (initializing_specie)
        {
          int status;

          /* Add current molecule component to link table */
          status = add_connection(specie_rule_indexes[current_specie_index],
                                  $1, current_object_index);
          if (FALSE == status)
          {
            yyerror_prev("Error adding link modifier");
            assert(0);
            exit(1);
          }
        }
        else
        {
          /* Add current molecule component rule to link table */
          status = add_connection(specie_rule_indexes[current_specie_index],
                                  $1, current_object_rule_index);
          if (FALSE == status)
          {
            yyerror_prev("Error adding link modifier");
            assert(0);
            exit(1);
          }
        }

        $$ = -1;
      }
  ;

molecule_component_modifier:
    '!' molecule_link_modifier
      {
        $$ = $2;
      }
  | object_modifier
      {
        $$ = $1;
      }
  ;

object_modifier:
    '~' name_or_number_or_wildcard
      {
        int status;
        color_t *state_color;

        state_color = real_dyncolors_color_create($2, COLOR_TYPE_STATE,
                                                  world->dynspecies,
                                                  world->dyncolors, world);
        if (NULL == state_color)
        {
            yyerror_prev("Error creating new color");
            assert(0);
            exit(1);
        }

        if (initializing_specie)
        {
          dynints_t *state_dynints;

          if (COLOR_TYPE_WILDCARD == state_color->type)
          {
            yyerror_prev("Not allowed to wildcard state in species.");
            exit(1);
          }

          /* Add state to current object's node */
          /* Store state list pointer for convenience. */
          state_dynints = current_specie_cgraph->
                          state_dynints_array.a[current_object_index];

          if (state_dynints->count > 0)

          if (SPECIES_BLOCK_TYPE == current_block)
          {
            // Only allow multipe states in molecule types block
            // or reaction rules.
            status = dynints_int_exists(state_color->value, state_dynints);
            {
              yyerror_prev("Not allowed to have more than one state per "
                           "component.");
              exit(1);
            }
          }

          status = dynints_int_add(state_color->value, state_dynints);
          if (FALSE == status)
          {
            yyerror_prev("Error adding object state");
            assert(0);
            exit(1);
          }
        }
        else
        {
          rule_t *object_rule;

          /* Add state to current object's rule */
          /* Store rule pointer for convenience. */
          object_rule = world->dynrules->rules.a + current_object_rule_index;

          status = dynints_int_exists(state_color->value,
                                      object_rule->state_dynints);
          if (TRUE == status)
          {
            yyerror_prev("Not allowed to list same state more than once.");
            assert(0);
            exit(1);
          }

          /* Make note of state. */
          status = dynints_int_add(state_color->value,
                                   object_rule->state_dynints);
          if (FALSE == status)
          {
            yyerror_prev("Error adding object state");
            assert(0);
            exit(1);
          }
        }
        $$ = -1;
      }
  | '@' name_or_number_or_wildcard
      {
        yyerror_prev("Partitions not implemented yet");
        assert(0);
        exit(1);

        $$ = -1;
      }
  | '=' expression_element
      {
        rule_t *object_rule;

        if (TRUE == initializing_specie)
        {
          yyerror_prev("Equals operator not yet implemented in "
                  "specie initialization yet");
          assert(0);
          exit(1);
        }
        else
        {
          /* Store rule pointer for convenience. */
          object_rule = world->dynrules->rules.a + current_object_rule_index;

          /* Make note of component multiple. */
          object_rule->count_flag = EQUAL;
          object_rule->count = $2;
        }

        /* Pass multiple up the parse tree. */
        $$ = $2;
      }
  | CMP_LT expression_element
      {
        rule_t *object_rule;

        if (TRUE == initializing_specie)
        {
          yyerror_prev("Comparison operators not allowed in specie initialization");
          assert(0);
          exit(1);
        }
        else
        {
          /* Store rule pointer for convenience. */
          object_rule = world->dynrules->rules.a + current_object_rule_index;

          /* Make note of component multiple restriction. */
          object_rule->count_flag = LESS_EQ;
          object_rule->count = $2 - 1;
        }

        $$ = -1;
      }
  | CMP_GT expression_element
      {
        rule_t *object_rule;

        if (TRUE == initializing_specie)
        {
          yyerror_prev("Comparison operators not allowed in specie initialization");
          assert(0);
          exit(1);
        }
        else
        {
          /* Store rule pointer for convenience. */
          object_rule = world->dynrules->rules.a + current_object_rule_index;

          /* Make note of component multiple restriction. */
          object_rule->count_flag = GREAT_EQ;
          object_rule->count = $2 + 1;
        }

        $$ = -1;
      }
  | CMP_LE expression_element
      {
        rule_t *object_rule;

        if (TRUE == initializing_specie)
        {
          yyerror_prev("Comparison operators not allowed in specie initialization");
          assert(0);
          exit(1);
        }
        else
        {
          /* Store rule pointer for convenience. */
          object_rule = world->dynrules->rules.a + current_object_rule_index;

          /* Make note of component multiple restriction. */
          object_rule->count_flag = LESS_EQ;
          object_rule->count = $2;
        }

        $$ = -1;
      }
  | CMP_GE expression_element
      {
        rule_t *object_rule;

        if (TRUE == initializing_specie)
        {
          yyerror_prev("Comparison operators not allowed in specie initialization");
          assert(0);
          exit(1);
        }
        else
        {
          /* Store rule pointer for convenience. */
          object_rule = world->dynrules->rules.a + current_object_rule_index;

          /* Make note of component multiple restriction. */
          object_rule->count_flag = GREAT_EQ;
          object_rule->count = $2;
        }

        $$ = -1;
      }
  ;

name_or_number:
    NAME
      { strcpy($$, $1); }
  | NUMBER
      { strcpy($$, $1); }
  ;

name_or_number_or_wildcard:
    NAME
      { strcpy($$, $1); }
  | NUMBER
      { strcpy($$, $1); }
  | '*'
      { strcpy($$, "*"); }
  ;

optional_label_and_or_number:
    /* Empty */
      {
        *$$ = '\0';
      }
  | NUMBER
      {
        strcpy($$, $1);
      }
  | QSTRING
      {
        strcpy($$, $1);
      }
  | LABEL
      {
        strcpy($$, $1);
      }
  | LABEL NUMBER
      {
        strcpy($$, $1);
      }
  ;

perl_hash:
    '{' perl_hash_item_list '}'
      {
        /* Finalize current hash. */
        ;
      }
  ;

perl_hash_item_list:
    /* Empty perl hash */
  | perl_hash_item
      {
        /* Create new hash structure and add current item to hash. */
        ;
      }
  | perl_hash_item_list ',' perl_hash_item
      {
        /* Add current item to existing hash. */
        ;
      }
  ;

perl_hash_item:
    NAME ',' perl_hash
      {
        /* Store current item. */
        yyerror_prev("Do not support hashes of hashes");
        assert(0);
        exit(1);
      }
  | NAME ',' perl_array
      {
        /* Store current item. */
        yyerror_prev("Do not support hashes of arrays");
        assert(0);
        exit(1);
      }
  | NAME ',' expression
      {
        /* Store current number. */
        bng_command_set_double($1, $3);
      }
  | NAME ',' QSTRING
      {
        /* Store current item. */
        bng_command_set_string($1, $3);
      }
  ;

perl_array_value:
    perl_hash
      {
        /* Store current item. */
        yyerror_prev("Do not support arrays of hashes");
        assert(0);
        exit(1);
      }
  | perl_array
      {
        /* Store current item. */
        yyerror_prev("Do not support arrays of arrays");
        assert(0);
        exit(1);
      }
  | expression
      {
        /* Store current number. */
        yyerror_prev("Do not support arrays of expressions");
        assert(0);
        exit(1);
      }
  | QSTRING
      {
        /* Store current item. */
        yyerror_prev("Do not support arrays of strings");
        assert(0);
        exit(1);
      }
  ;

perl_array:
    '[' perl_array_item_list ']'
      {
        /* Finalize current array. */
        ;
      }
  ;

perl_array_item_list:
    /* Empty perl array */
  | perl_array_value
      { 
        ; 
      }
  | perl_array_item_list ',' perl_array_value
      { 
        ; 
      }
  ;

expression: 
    expression_element { $$ = $1; }
  | expression '+' expression { $$ = $1 + $3; }
  | expression '-' expression { $$ = $1 - $3; }
  | expression '*' expression { $$ = $1 * $3; }
  | expression '/' expression
      {
        if ($3 == 0.0)
        {
          yyerror_prev("divide by zero");
          assert(0);
          exit(1);
        }
        else
        {
          $$ = $1 / $3;
        }
      }
  | '-' expression %prec UMINUS { $$ = -$2; }
  | expression '^' expression { $$ = pow($1, $3); }
  | '(' expression ')' { $$ = $2; }
  ;

expression_element:
    NUMBER
      {
        char *endptr;

        $$ = strtod($1, &endptr);
        if (endptr == $1)
        {
          yyerror_prev("Error occured when trying to convert string to number");
          assert(0);
          exit(1);
        }
      }
  | symbol
      {
        size_t index;
        double value;

        if (0 == strncasecmp($1, SYSTEM_PREFIX, sizeof(SYSTEM_PREFIX) - 1))
        {
          char *system_symbol;

          system_symbol = $1;

          // Skip over prefix.
          system_symbol += sizeof(SYSTEM_PREFIX) - 1;

          index = ds_system_get_value(system_symbol, &value);
          if (0 == index)
          {
            yyerror("system symbol not found");
            assert(0);
            exit(1);
          }
        }
        else
        {
          index = dynsymbols_get_value($1, &value, dynsymbols);
          if (0 == index)
          {
            yyerror("symbol not found");
            assert(0);
            exit(1);
          }
        }

        $$ = value;
      }
  ;


%%

