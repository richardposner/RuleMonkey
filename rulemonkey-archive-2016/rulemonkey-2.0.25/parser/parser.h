/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     LABEL = 258,
     QSTRING = 259,
     NAME = 260,
     NUMBER = 261,
     QUESTION_MARK = 262,
     MOLECULES = 263,
     SPECIES = 264,
     BEGIN_PARAMETERS = 265,
     END_PARAMETERS = 266,
     BEGIN_SPECIES = 267,
     END_SPECIES = 268,
     BEGIN_REACTION_RULES = 269,
     END_REACTION_RULES = 270,
     BEGIN_OBSERVABLES = 271,
     END_OBSERVABLES = 272,
     BEGIN_MOLECULE_TYPES = 273,
     END_MOLECULE_TYPES = 274,
     BEGIN_ACTIONS = 275,
     END_ACTIONS = 276,
     UNIDIRECTIONAL_RULE = 277,
     BIDIRECTIONAL_RULE = 278,
     GENERATE_NETWORK = 279,
     SIMULATE_ODE = 280,
     SIMULATE_SSA = 281,
     TO_SBML_FILE = 282,
     SETUP_PS = 283,
     SIMULATE_PS = 284,
     SIMULATE_RS = 285,
     STATE_OUTPUT_PS = 286,
     CMP_LE = 287,
     CMP_GE = 288,
     CMP_LT = 289,
     CMP_GT = 290,
     LOGICAL_AND = 291,
     LOGICAL_OR = 292,
     LOGICAL_XOR = 293,
     UMINUS = 294
   };
#endif
/* Tokens.  */
#define LABEL 258
#define QSTRING 259
#define NAME 260
#define NUMBER 261
#define QUESTION_MARK 262
#define MOLECULES 263
#define SPECIES 264
#define BEGIN_PARAMETERS 265
#define END_PARAMETERS 266
#define BEGIN_SPECIES 267
#define END_SPECIES 268
#define BEGIN_REACTION_RULES 269
#define END_REACTION_RULES 270
#define BEGIN_OBSERVABLES 271
#define END_OBSERVABLES 272
#define BEGIN_MOLECULE_TYPES 273
#define END_MOLECULE_TYPES 274
#define BEGIN_ACTIONS 275
#define END_ACTIONS 276
#define UNIDIRECTIONAL_RULE 277
#define BIDIRECTIONAL_RULE 278
#define GENERATE_NETWORK 279
#define SIMULATE_ODE 280
#define SIMULATE_SSA 281
#define TO_SBML_FILE 282
#define SETUP_PS 283
#define SIMULATE_PS 284
#define SIMULATE_RS 285
#define STATE_OUTPUT_PS 286
#define CMP_LE 287
#define CMP_GE 288
#define CMP_LT 289
#define CMP_GT 290
#define LOGICAL_AND 291
#define LOGICAL_OR 292
#define LOGICAL_XOR 293
#define UMINUS 294




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 34 "parser.y"
{
  int ival;
  double dval;
  char str[1024];
}
/* Line 1529 of yacc.c.  */
#line 133 "parser.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

