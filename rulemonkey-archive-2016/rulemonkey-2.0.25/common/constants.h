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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "../nauty22/nauty.h"

#define DEFAULT_DEBUG_LEVEL 0

#define AVOGADRO 6.0221415e23

#define SYSTEM_PREFIX "system_"

typedef enum reactant_type_enum reactant_type_t;
enum reactant_type_enum
{
  REACTANT_TYPE_A,
  REACTANT_TYPE_B,
  REACTANT_TYPE_AA,
  REACTANT_TYPE_AB
};

//#define EXERCISE_ALLOC

#define MAX_BOOLEAN_STR_LENGTH (sizeof("FALSE") + 1)

#define MAX_REACTANTS 2

// Maximum number of reactants + maximum number of products.
// 2 + 1 = 3
// Make it larger to account for multiple species in a given observable.
#define MAX_SPECIE_RULE_INDEX 20

#ifdef EXERCISE_ALLOC
#define INITIAL_GRAPH_SIZE 1
#define INITIAL_DYNCOMBINATIONS 1
#define INITIAL_DYNSPECIES 1
#define INITIAL_DYNSPECIE_CONTAINERS 1
#define INITIAL_DYNRULES 1
#define INITIAL_DYNCOLORS 1
#define INITIAL_DYNREACTIONS 1
#define INITIAL_DYNREACTION_POINTERS 1
#define INITIAL_DYNREACTION_CONTAINERS 1
#define INITIAL_MOLECULE_COMPONENTS 1
#define INITIAL_DYNSYMBOLS 1
#define INITIAL_DYNMOLECULES 1
#define INITIAL_CGRAPH_SIZE 1
#define INITIAL_DYNCONNECTIONS 1
#define INITIAL_DYNINTS 1
#define INITIAL_DYNMAPS 1
#define INITIAL_RULE_AGGREGATES 1
#define INITIAL_REACTION_CONTAINERS 1
#define INITIAL_DYNPS_SPECIES 1
#define INITIAL_DYNPS_SPECIES 1
#define INITIAL_DYNREPORT_DEFINITION_CONTAINERS 1
#define INITIAL_FSA_BLOCK_COUNT 1
#define INITIAL_SPECIES_FSA_BLOCK_COUNT 1
#define INITIAL_STATES 1
#define INITIAL_COMP_PERM_COUNT 1
#else
#define INITIAL_GRAPH_SIZE 16
#define INITIAL_DYNCOMBINATIONS 2
#define INITIAL_DYNSPECIES 50000
#define INITIAL_DYNSPECIE_CONTAINERS 8
#define INITIAL_DYNRULES 16
#define INITIAL_DYNCOLORS 16
#define INITIAL_DYNREACTIONS 8192
#define INITIAL_DYNREACTION_POINTERS 1024
#define INITIAL_DYNREACTION_CONTAINERS 16
#define INITIAL_MOLECULE_COMPONENTS 4
#define INITIAL_DYNSYMBOLS 16
#define INITIAL_DYNMOLECULES 16
#define INITIAL_CGRAPH_SIZE 16
#define INITIAL_DYNCONNECTIONS 4
#define INITIAL_DYNINTS 2
#define INITIAL_DYNMAPS 8
#define INITIAL_RULE_AGGREGATES 4
#define INITIAL_REACTION_CONTAINERS 8
#define INITIAL_DYNPS_SPECIES INITIAL_DYNSPECIES
#define INITIAL_DYNREPORT_DEFINITION_CONTAINERS 8
#define INITIAL_FSA_BLOCK_COUNT 4
#define INITIAL_STATES 4
#define INITIAL_COMP_PERM_COUNT 8
#endif

// To prevent warnings from nauty include files.
#ifdef SETWORD_INT
#define SUPPRESS_WARNINGS printf("%i%i%i", *bit, *bytecount, *leftbit)
#endif

#ifdef SETWORD_LONG
#define SUPPRESS_WARNINGS printf("%ld%i%i", *bit, *bytecount, *leftbit)
#endif

#ifdef SETWORD_LONGLONG
#define SUPPRESS_WARNINGS printf("%lld%i%i", *bit, *bytecount, *leftbit)
#endif


#endif
