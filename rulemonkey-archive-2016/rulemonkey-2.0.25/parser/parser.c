/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



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




/* Copy the first part of user declarations.  */
#line 2 "parser.y"

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



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 34 "parser.y"
{
  int ival;
  double dval;
  char str[1024];
}
/* Line 193 of yacc.c.  */
#line 212 "parser.c"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 225 "parser.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  32
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   226

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  58
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  79
/* YYNRULES -- Number of rules.  */
#define YYNRULES  160
/* YYNRULES -- Number of states.  */
#define YYNSTATES  231

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   294

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    51,     2,     2,     2,     2,     2,     2,
      45,    46,    41,    40,    49,    39,    50,    42,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    47,
       2,    48,     2,     2,    53,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    56,     2,    57,    44,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    54,     2,    55,    52,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    43
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     5,     7,    10,    12,    14,    16,    18,
      20,    22,    24,    30,    31,    33,    35,    37,    39,    41,
      43,    45,    47,    49,    52,    56,    59,    63,    66,    70,
      73,    77,    80,    84,    87,    91,    93,    95,    97,    99,
     101,   103,   105,   107,   109,   111,   113,   115,   117,   120,
     122,   125,   127,   130,   132,   135,   137,   140,   142,   145,
     148,   151,   155,   157,   159,   163,   165,   167,   169,   175,
     183,   185,   187,   189,   192,   194,   196,   200,   204,   209,
     214,   218,   220,   224,   228,   232,   234,   236,   238,   240,
     245,   247,   251,   253,   257,   260,   263,   266,   268,   270,
     274,   275,   277,   281,   284,   286,   287,   289,   291,   294,
     295,   297,   299,   302,   304,   306,   308,   311,   313,   316,
     319,   322,   325,   328,   331,   334,   336,   338,   340,   342,
     344,   345,   347,   349,   351,   354,   358,   359,   361,   365,
     369,   373,   377,   381,   383,   385,   387,   389,   393,   394,
     396,   400,   402,   406,   410,   414,   418,   421,   425,   429,
     431
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
      59,     0,    -1,    60,    -1,    61,    -1,    60,    61,    -1,
      62,    -1,    65,    -1,    66,    -1,    67,    -1,    68,    -1,
      69,    -1,    70,    -1,    64,    45,    63,    46,    47,    -1,
      -1,   129,    -1,    24,    -1,    25,    -1,    26,    -1,    31,
      -1,    28,    -1,    29,    -1,    30,    -1,    27,    -1,    71,
      72,    -1,    71,    83,    72,    -1,    73,    74,    -1,    73,
      84,    74,    -1,    75,    76,    -1,    75,    85,    76,    -1,
      77,    78,    -1,    77,    86,    78,    -1,    79,    80,    -1,
      79,    87,    80,    -1,    81,    82,    -1,    81,    88,    82,
      -1,    10,    -1,    11,    -1,    12,    -1,    13,    -1,    14,
      -1,    15,    -1,    16,    -1,    17,    -1,    18,    -1,    19,
      -1,    20,    -1,    21,    -1,    89,    -1,    83,    89,    -1,
      93,    -1,    84,    93,    -1,    95,    -1,    85,    95,    -1,
     102,    -1,    86,   102,    -1,    94,    -1,    87,    94,    -1,
      92,    -1,    88,    92,    -1,   128,    90,    -1,    91,   135,
      -1,    91,    48,   135,    -1,     5,    -1,    62,    -1,   128,
     105,   135,    -1,   112,    -1,    96,    -1,    97,    -1,   128,
      98,    22,    99,   135,    -1,   128,    98,    23,    99,   135,
      49,   135,    -1,   101,    -1,   101,    -1,   103,    -1,   100,
     103,    -1,    41,    -1,   103,    -1,   101,    40,    41,    -1,
     101,    40,   103,    -1,   128,     8,     5,   100,    -1,   128,
       9,     5,   100,    -1,    45,   103,    46,    -1,   104,    -1,
     103,    36,   104,    -1,   103,    37,   104,    -1,   103,    38,
     104,    -1,   106,    -1,   107,    -1,   109,    -1,   108,    -1,
      45,   108,    46,   121,    -1,   112,    -1,   108,    50,   112,
      -1,   110,    -1,   109,    50,   110,    -1,   111,   121,    -1,
     113,   115,    -1,   114,   115,    -1,     5,    -1,     5,    -1,
      45,   116,    46,    -1,    -1,   117,    -1,   116,    49,   117,
      -1,   118,   119,    -1,     5,    -1,    -1,   120,    -1,   124,
      -1,   120,   124,    -1,    -1,   122,    -1,   125,    -1,   122,
     125,    -1,    40,    -1,     7,    -1,   126,    -1,    51,   123,
      -1,   125,    -1,    52,   127,    -1,    53,   127,    -1,    48,
     136,    -1,    34,   136,    -1,    35,   136,    -1,    32,   136,
      -1,    33,   136,    -1,     5,    -1,     6,    -1,     5,    -1,
       6,    -1,    41,    -1,    -1,     6,    -1,     4,    -1,     3,
      -1,     3,     6,    -1,    54,   130,    55,    -1,    -1,   131,
      -1,   130,    49,   131,    -1,     5,    49,   129,    -1,     5,
      49,   133,    -1,     5,    49,   135,    -1,     5,    49,     4,
      -1,   129,    -1,   133,    -1,   135,    -1,     4,    -1,    56,
     134,    57,    -1,    -1,   132,    -1,   134,    49,   132,    -1,
     136,    -1,   135,    40,   135,    -1,   135,    39,   135,    -1,
     135,    41,   135,    -1,   135,    42,   135,    -1,    39,   135,
      -1,   135,    44,   135,    -1,    45,   135,    46,    -1,     6,
      -1,    91,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    93,    93,    98,   100,   105,   107,   109,   111,   113,
     115,   117,   122,   139,   141,   146,   151,   156,   161,   168,
     175,   190,   197,   205,   207,   212,   214,   219,   221,   226,
     228,   233,   235,   240,   242,   247,   252,   257,   262,   267,
     272,   277,   282,   287,   294,   302,   307,   312,   314,   319,
     321,   326,   328,   333,   335,   340,   342,   347,   349,   354,
     360,   391,   425,   429,   434,   499,   512,   517,   525,   538,
     552,   566,   579,   586,   596,   603,   611,   619,   627,   642,
     661,   666,   671,   678,   686,   697,   708,   719,   736,   740,
     751,   756,   766,   771,   781,   797,   808,   819,   894,   926,
     938,   940,   942,   947,  1003,  1071,  1074,  1081,  1092,  1112,
    1115,  1122,  1133,  1152,  1174,  1196,  1233,  1237,  1244,  1325,
    1333,  1357,  1379,  1401,  1423,  1448,  1450,  1455,  1457,  1459,
    1465,  1468,  1472,  1476,  1480,  1487,  1494,  1496,  1501,  1509,
    1516,  1523,  1528,  1536,  1543,  1550,  1557,  1567,  1574,  1576,
    1580,  1587,  1588,  1589,  1590,  1591,  1604,  1605,  1606,  1610,
    1622
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "LABEL", "QSTRING", "NAME", "NUMBER",
  "QUESTION_MARK", "MOLECULES", "SPECIES", "BEGIN_PARAMETERS",
  "END_PARAMETERS", "BEGIN_SPECIES", "END_SPECIES", "BEGIN_REACTION_RULES",
  "END_REACTION_RULES", "BEGIN_OBSERVABLES", "END_OBSERVABLES",
  "BEGIN_MOLECULE_TYPES", "END_MOLECULE_TYPES", "BEGIN_ACTIONS",
  "END_ACTIONS", "UNIDIRECTIONAL_RULE", "BIDIRECTIONAL_RULE",
  "GENERATE_NETWORK", "SIMULATE_ODE", "SIMULATE_SSA", "TO_SBML_FILE",
  "SETUP_PS", "SIMULATE_PS", "SIMULATE_RS", "STATE_OUTPUT_PS", "CMP_LE",
  "CMP_GE", "CMP_LT", "CMP_GT", "LOGICAL_AND", "LOGICAL_OR", "LOGICAL_XOR",
  "'-'", "'+'", "'*'", "'/'", "UMINUS", "'^'", "'('", "')'", "';'", "'='",
  "','", "'.'", "'!'", "'~'", "'@'", "'{'", "'}'", "'['", "']'", "$accept",
  "bngl_file", "bng_item_list", "bng_item", "bng_command",
  "optional_perl_hash", "bng_command_name", "parameter_block",
  "specie_block", "reaction_rule_block", "observable_block",
  "molecule_types_block", "action_block", "begin_parameters",
  "end_parameters", "begin_species", "end_species", "begin_reaction_rules",
  "end_reaction_rules", "begin_observables", "end_observables",
  "begin_molecule_types", "end_molecule_types", "begin_actions",
  "end_actions", "parameter_block_list", "specie_block_list",
  "reaction_rule_block_list", "observable_block_list",
  "molecule_type_block_list", "action_block_list",
  "parameter_block_list_item", "parameter_declaration", "symbol",
  "action_block_list_item", "specie_block_list_item",
  "molecule_type_block_list_item", "reaction_rule_block_list_item",
  "reaction_rule", "reaction_rule_type", "reaction_rule_lhs",
  "reaction_rule_rhs", "specie_rule_or_list", "specie_rule_plus_list",
  "observable_block_list_item", "specie_rule", "specie",
  "specie_declaration", "molecule_aggregate_modified",
  "molecule_aggregate_modified_declaration",
  "molecule_aggregate_declaration", "molecule_aggregate",
  "molecule_modified", "molecule", "molecule_declaration", "molecule_name",
  "molecule_name_declaration", "molecule_components",
  "molecule_component_list", "molecule_component",
  "molecule_component_name", "optional_molecule_component_modifier_list",
  "molecule_component_modifier_list", "optional_object_modifier_list",
  "object_modifier_list", "molecule_link_modifier",
  "molecule_component_modifier", "object_modifier", "name_or_number",
  "name_or_number_or_wildcard", "optional_label_and_or_number",
  "perl_hash", "perl_hash_item_list", "perl_hash_item", "perl_array_value",
  "perl_array", "perl_array_item_list", "expression", "expression_element", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,    45,
      43,    42,    47,   294,    94,    40,    41,    59,    61,    44,
      46,    33,   126,    64,   123,   125,    91,    93
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    58,    59,    60,    60,    61,    61,    61,    61,    61,
      61,    61,    62,    63,    63,    64,    64,    64,    64,    64,
      64,    64,    64,    65,    65,    66,    66,    67,    67,    68,
      68,    69,    69,    70,    70,    71,    72,    73,    74,    75,
      76,    77,    78,    79,    80,    81,    82,    83,    83,    84,
      84,    85,    85,    86,    86,    87,    87,    88,    88,    89,
      90,    90,    91,    92,    93,    94,    95,    96,    97,    97,
      98,    99,   100,   100,   101,   101,   101,   101,   102,   102,
     103,   103,   103,   103,   103,   104,   105,   106,   107,   107,
     108,   108,   109,   109,   110,   111,   112,   113,   114,   115,
     116,   116,   116,   117,   118,   119,   119,   120,   120,   121,
     121,   122,   122,   123,   123,   123,   124,   124,   125,   125,
     125,   125,   125,   125,   125,   126,   126,   127,   127,   127,
     128,   128,   128,   128,   128,   129,   130,   130,   130,   131,
     131,   131,   131,   132,   132,   132,   132,   133,   134,   134,
     134,   135,   135,   135,   135,   135,   135,   135,   135,   136,
     136
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     2,     1,     1,     1,     1,     1,
       1,     1,     5,     0,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     2,     3,     2,     3,     2,     3,     2,
       3,     2,     3,     2,     3,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     2,     1,
       2,     1,     2,     1,     2,     1,     2,     1,     2,     2,
       2,     3,     1,     1,     3,     1,     1,     1,     5,     7,
       1,     1,     1,     2,     1,     1,     3,     3,     4,     4,
       3,     1,     3,     3,     3,     1,     1,     1,     1,     4,
       1,     3,     1,     3,     2,     2,     2,     1,     1,     3,
       0,     1,     3,     2,     1,     0,     1,     1,     2,     0,
       1,     1,     2,     1,     1,     1,     2,     1,     2,     2,
       2,     2,     2,     2,     2,     1,     1,     1,     1,     1,
       0,     1,     1,     1,     2,     3,     0,     1,     3,     3,
       3,     3,     3,     1,     1,     1,     1,     3,     0,     1,
       3,     1,     3,     3,     3,     3,     2,     3,     3,     1,
       1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,    35,    37,    39,    41,    43,    45,    15,    16,    17,
      22,    19,    20,    21,    18,     0,     2,     3,     5,     0,
       6,     7,     8,     9,    10,    11,   130,   130,   130,   130,
       0,     0,     1,     4,    13,   133,   132,   131,    36,    23,
     130,    47,     0,    38,    25,   130,    49,     0,    40,    27,
     130,    51,    66,    67,     0,    42,    29,   130,    53,     0,
      98,    44,    31,     0,    55,    65,     0,    46,    63,    33,
       0,    57,   136,     0,    14,   134,    24,    48,    62,    59,
       0,    26,    50,     0,     0,    86,    88,    90,    28,    52,
      97,    74,     0,     0,    70,    75,    81,    85,    87,    92,
     109,     0,    30,    54,     0,     0,    32,    56,   100,    96,
      34,    58,     0,     0,   137,     0,   159,     0,     0,     0,
     160,    60,   151,     0,    64,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    94,   110,   111,    95,     0,     0,   104,     0,   101,
     105,     0,     0,   135,    12,   156,     0,    61,     0,     0,
       0,     0,     0,   109,    91,    80,     0,    71,     0,    76,
      77,    82,    83,    84,    93,   123,   124,   121,   122,   120,
     127,   128,   129,   118,   119,   112,    78,    72,    79,    99,
       0,     0,   103,   106,   107,   117,   142,   148,   139,   140,
     141,   138,   158,   153,   152,   154,   155,   157,    89,    68,
       0,    73,   102,   125,   126,   114,   113,   116,   115,   108,
     146,   143,   149,   144,     0,   145,     0,     0,   147,    69,
     150
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    15,    16,    17,    18,    73,    19,    20,    21,    22,
      23,    24,    25,    26,    39,    27,    44,    28,    49,    29,
      56,    30,    62,    31,    69,    40,    45,    50,    57,    63,
      70,    41,    79,   120,    71,    46,    64,    51,    52,    53,
      93,   166,   186,   167,    58,    95,    96,    84,    97,    85,
      86,    98,    99,   100,    65,   101,    66,   109,   148,   149,
     150,   192,   193,   141,   142,   217,   194,   143,   218,   183,
      42,   221,   113,   114,   222,   223,   224,   225,   122
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -129
static const yytype_int16 yypact[] =
{
     110,  -129,  -129,  -129,  -129,  -129,  -129,  -129,  -129,  -129,
    -129,  -129,  -129,  -129,  -129,    50,   110,  -129,  -129,   -17,
    -129,  -129,  -129,  -129,  -129,  -129,   173,    97,    21,    81,
      16,   130,  -129,  -129,     0,    86,  -129,  -129,  -129,  -129,
     173,  -129,    99,  -129,  -129,    97,  -129,    -3,  -129,  -129,
      21,  -129,  -129,  -129,    10,  -129,  -129,    81,  -129,    58,
    -129,  -129,  -129,    16,  -129,  -129,    76,  -129,  -129,  -129,
     130,  -129,   106,    60,  -129,  -129,  -129,  -129,  -129,  -129,
      14,  -129,  -129,   120,    25,  -129,    77,  -129,  -129,  -129,
    -129,  -129,     2,    54,    89,    95,  -129,  -129,    94,  -129,
      61,    76,  -129,  -129,   147,   148,  -129,  -129,   157,  -129,
    -129,  -129,   116,   -26,  -129,   117,  -129,    25,    25,    25,
    -129,   154,  -129,    -2,   154,   120,   -20,    10,    10,    27,
     164,   164,   164,   164,   111,   111,   111,   111,   111,    28,
      28,  -129,    61,  -129,  -129,     2,     2,  -129,    56,  -129,
     115,     4,   106,  -129,  -129,   134,   146,   154,    25,    25,
      25,    25,    25,    61,  -129,  -129,    25,    89,    25,  -129,
      95,  -129,  -129,  -129,  -129,  -129,  -129,  -129,  -129,  -129,
    -129,  -129,  -129,  -129,  -129,  -129,     2,    95,     2,  -129,
     157,    68,  -129,   115,  -129,  -129,  -129,     7,  -129,  -129,
     154,  -129,  -129,   101,   101,   134,   134,   134,  -129,   154,
     131,    95,  -129,  -129,  -129,  -129,  -129,  -129,  -129,  -129,
    -129,  -129,  -129,  -129,    66,   154,    25,     7,  -129,   154,
    -129
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -129,  -129,  -129,   158,   -25,  -129,  -129,  -129,  -129,  -129,
    -129,  -129,  -129,  -129,   149,  -129,   152,  -129,   141,  -129,
     150,  -129,   142,  -129,   136,  -129,  -129,  -129,  -129,  -129,
    -129,   168,  -129,   162,   139,   165,   151,   161,  -129,  -129,
    -129,    84,    67,   163,   159,   -89,    51,  -129,  -129,  -129,
     132,  -129,    85,  -129,   -42,  -129,  -129,   118,  -129,    30,
    -129,  -129,  -129,    59,  -129,  -129,    31,  -128,  -129,    83,
      62,   -33,  -129,    69,    -1,    74,  -129,   -80,    65
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
     121,    74,    60,   126,   124,    87,    68,    90,   196,    78,
     116,   220,    78,   116,   185,    90,   130,   131,   132,    78,
     116,    60,   195,   152,    35,    36,   165,    37,    34,   153,
      78,   116,    90,   180,   181,    61,    48,   155,   156,   157,
     170,    87,    83,   117,   163,    68,   117,    92,   125,   118,
      32,    91,   118,   117,    72,    92,   187,   187,    72,   118,
     197,    72,   119,   197,   117,   195,   104,   105,   169,   182,
     118,   200,    92,   213,   214,   215,   127,   128,   203,   204,
     205,   206,   207,   164,    35,    36,   209,    37,   210,    47,
      54,    59,    75,   134,   135,   136,   137,   211,    55,   211,
      35,    36,   189,    37,    78,   190,   115,    47,   216,   138,
      43,   112,    54,   139,   140,   227,    78,   116,   198,    59,
       1,   108,     2,   228,     3,    60,     4,   125,     5,   129,
       6,   130,   131,   132,     7,     8,     9,    10,    11,    12,
      13,    14,   160,   161,   133,   162,   229,   134,   135,   136,
     137,    67,   145,   146,     7,     8,     9,    10,    11,    12,
      13,    14,   147,   138,   154,   151,   191,   139,   140,    90,
     158,   159,   160,   161,    33,   162,    35,    36,   162,    37,
     226,   171,   172,   173,    38,   158,   159,   160,   161,    76,
     162,    88,   202,   158,   159,   160,   161,    81,   162,   175,
     176,   177,   178,   179,    80,   106,   110,   102,    77,   111,
      82,    89,   168,   188,   107,   123,   103,    94,   174,   144,
     212,   201,   208,   184,   219,   199,   230
};

static const yytype_uint8 yycheck[] =
{
      80,    34,     5,    92,    84,    47,    31,     5,     4,     5,
       6,     4,     5,     6,   142,     5,    36,    37,    38,     5,
       6,     5,   150,    49,     3,     4,    46,     6,    45,    55,
       5,     6,     5,     5,     6,    19,    15,   117,   118,   119,
     129,    83,    45,    39,    46,    70,    39,    45,    50,    45,
       0,    41,    45,    39,    54,    45,   145,   146,    54,    45,
      56,    54,    48,    56,    39,   193,     8,     9,    41,    41,
      45,   151,    45,     5,     6,     7,    22,    23,   158,   159,
     160,   161,   162,   125,     3,     4,   166,     6,   168,    27,
      28,    29,     6,    32,    33,    34,    35,   186,    17,   188,
       3,     4,    46,     6,     5,    49,    46,    45,    40,    48,
      13,     5,    50,    52,    53,    49,     5,     6,   151,    57,
      10,    45,    12,    57,    14,     5,    16,    50,    18,    40,
      20,    36,    37,    38,    24,    25,    26,    27,    28,    29,
      30,    31,    41,    42,    50,    44,   226,    32,    33,    34,
      35,    21,     5,     5,    24,    25,    26,    27,    28,    29,
      30,    31,     5,    48,    47,    49,    51,    52,    53,     5,
      39,    40,    41,    42,    16,    44,     3,     4,    44,     6,
      49,   130,   131,   132,    11,    39,    40,    41,    42,    40,
      44,    50,    46,    39,    40,    41,    42,    45,    44,   134,
     135,   136,   137,   138,    42,    63,    70,    57,    40,    70,
      45,    50,   128,   146,    63,    83,    57,    54,   133,   101,
     190,   152,   163,   140,   193,   151,   227
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    10,    12,    14,    16,    18,    20,    24,    25,    26,
      27,    28,    29,    30,    31,    59,    60,    61,    62,    64,
      65,    66,    67,    68,    69,    70,    71,    73,    75,    77,
      79,    81,     0,    61,    45,     3,     4,     6,    11,    72,
      83,    89,   128,    13,    74,    84,    93,   128,    15,    76,
      85,    95,    96,    97,   128,    17,    78,    86,   102,   128,
       5,    19,    80,    87,    94,   112,   114,    21,    62,    82,
      88,    92,    54,    63,   129,     6,    72,    89,     5,    90,
      91,    74,    93,    45,   105,   107,   108,   112,    76,    95,
       5,    41,    45,    98,   101,   103,   104,   106,   109,   110,
     111,   113,    78,   102,     8,     9,    80,    94,    45,   115,
      82,    92,     5,   130,   131,    46,     6,    39,    45,    48,
      91,   135,   136,   108,   135,    50,   103,    22,    23,    40,
      36,    37,    38,    50,    32,    33,    34,    35,    48,    52,
      53,   121,   122,   125,   115,     5,     5,     5,   116,   117,
     118,    49,    49,    55,    47,   135,   135,   135,    39,    40,
      41,    42,    44,    46,   112,    46,    99,   101,    99,    41,
     103,   104,   104,   104,   110,   136,   136,   136,   136,   136,
       5,     6,    41,   127,   127,   125,   100,   103,   100,    46,
      49,    51,   119,   120,   124,   125,     4,    56,   129,   133,
     135,   131,    46,   135,   135,   135,   135,   135,   121,   135,
     135,   103,   117,     5,     6,     7,    40,   123,   126,   124,
       4,   129,   132,   133,   134,   135,    49,    49,    57,   135,
     132
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 94 "parser.y"
    { ; }
    break;

  case 3:
#line 99 "parser.y"
    { ; }
    break;

  case 4:
#line 101 "parser.y"
    { ; }
    break;

  case 5:
#line 106 "parser.y"
    { ; }
    break;

  case 6:
#line 108 "parser.y"
    { ; }
    break;

  case 7:
#line 110 "parser.y"
    { ; }
    break;

  case 8:
#line 112 "parser.y"
    { ; }
    break;

  case 9:
#line 114 "parser.y"
    { ; }
    break;

  case 10:
#line 116 "parser.y"
    { ; }
    break;

  case 11:
#line 118 "parser.y"
    { ; }
    break;

  case 12:
#line 123 "parser.y"
    {
        boolean status;

        if (*(yyvsp[(1) - (5)].str) != '\0')
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
    break;

  case 14:
#line 142 "parser.y"
    { ; }
    break;

  case 15:
#line 147 "parser.y"
    {
        yyerror_prev("generate_network() not implemented.");
        *(yyval.str) = '\0';
      }
    break;

  case 16:
#line 152 "parser.y"
    {
        yyerror_prev("simulate_ode() not implemented.");
        *(yyval.str) = '\0';
      }
    break;

  case 17:
#line 157 "parser.y"
    {
        yyerror_prev("simulate_ssa() not implemented.");
        *(yyval.str) = '\0';
      }
    break;

  case 18:
#line 162 "parser.y"
    {
        /* Set callback. */
        bng_command_ptr = state_output_ps;
        bng_command_set_default();
        strcpy((yyval.str), "ps_state_output");
      }
    break;

  case 19:
#line 169 "parser.y"
    {
        /* Set callback. */
        bng_command_ptr = setup_ps;
        bng_command_set_default();
        strcpy((yyval.str), "setup_ps");
      }
    break;

  case 20:
#line 176 "parser.y"
    {
// simulate_ps does not currently work, rule checking has infinite loop.
#ifdef ENABLE_SIMULATE_PS
        // simulate_ps does not currently work, rule checking has infinite loop.
        /* Set callback. */
        bng_command_ptr = simulate_ps;
        bng_command_set_default();
        strcpy((yyval.str), "simulate_ps");
#else
        // simulate_ps does not currently work, rule checking has infinite loop.
        yyerror_prev("simulate_ps() not implemented.");
        *(yyval.str) = '\0';
#endif
      }
    break;

  case 21:
#line 191 "parser.y"
    {
        /* Set callback. */
        bng_command_ptr = simulate_rs;
        bng_command_set_default();
        strcpy((yyval.str), "simulate_rs");
      }
    break;

  case 22:
#line 198 "parser.y"
    {
        yyerror_prev("to_sbml_file() not implemented.");
        *(yyval.str) = '\0';
      }
    break;

  case 23:
#line 206 "parser.y"
    { ; }
    break;

  case 24:
#line 208 "parser.y"
    { ; }
    break;

  case 25:
#line 213 "parser.y"
    { ; }
    break;

  case 26:
#line 215 "parser.y"
    { ; }
    break;

  case 27:
#line 220 "parser.y"
    { ; }
    break;

  case 28:
#line 222 "parser.y"
    { ; }
    break;

  case 29:
#line 227 "parser.y"
    { ; }
    break;

  case 30:
#line 229 "parser.y"
    { ; }
    break;

  case 31:
#line 234 "parser.y"
    { ; }
    break;

  case 32:
#line 236 "parser.y"
    { ; }
    break;

  case 33:
#line 241 "parser.y"
    { ; }
    break;

  case 34:
#line 243 "parser.y"
    { ; }
    break;

  case 35:
#line 248 "parser.y"
    { current_block = PARAMETERS_BLOCK_TYPE; }
    break;

  case 36:
#line 253 "parser.y"
    { current_block = NO_BLOCK_TYPE; }
    break;

  case 37:
#line 258 "parser.y"
    { current_block = SPECIES_BLOCK_TYPE; }
    break;

  case 38:
#line 263 "parser.y"
    { current_block = NO_BLOCK_TYPE; }
    break;

  case 39:
#line 268 "parser.y"
    { current_block = REACTION_RULES_BLOCK_TYPE; }
    break;

  case 40:
#line 273 "parser.y"
    { current_block = NO_BLOCK_TYPE; }
    break;

  case 41:
#line 278 "parser.y"
    { current_block = OBSERVABLES_BLOCK_TYPE; }
    break;

  case 42:
#line 283 "parser.y"
    { current_block = NO_BLOCK_TYPE; }
    break;

  case 43:
#line 288 "parser.y"
    {
        current_block = MOLECULE_TYPES_BLOCK_TYPE;
      }
    break;

  case 44:
#line 295 "parser.y"
    {
        current_block = NO_BLOCK_TYPE;
        molecule_types_present = 1;
      }
    break;

  case 45:
#line 303 "parser.y"
    { current_block = ACTIONS_BLOCK_TYPE; }
    break;

  case 46:
#line 308 "parser.y"
    { current_block = NO_BLOCK_TYPE; }
    break;

  case 47:
#line 313 "parser.y"
    { ; }
    break;

  case 48:
#line 315 "parser.y"
    { ; }
    break;

  case 49:
#line 320 "parser.y"
    { ; }
    break;

  case 50:
#line 322 "parser.y"
    { ; }
    break;

  case 51:
#line 327 "parser.y"
    { ; }
    break;

  case 52:
#line 329 "parser.y"
    { ; }
    break;

  case 53:
#line 334 "parser.y"
    { ; }
    break;

  case 54:
#line 336 "parser.y"
    { ; }
    break;

  case 55:
#line 341 "parser.y"
    { ; }
    break;

  case 56:
#line 343 "parser.y"
    { ; }
    break;

  case 57:
#line 348 "parser.y"
    { ; }
    break;

  case 58:
#line 350 "parser.y"
    { ; }
    break;

  case 59:
#line 355 "parser.y"
    {
      }
    break;

  case 60:
#line 361 "parser.y"
    {
        int status;

        if (0 == strncasecmp((yyvsp[(1) - (2)].str), SYSTEM_PREFIX, sizeof(SYSTEM_PREFIX) - 1))
        {
          char *system_symbol;

          system_symbol = (yyvsp[(1) - (2)].str);

          // Skip over prefix.
          system_symbol += sizeof(SYSTEM_PREFIX) - 1;

          status = ds_system_set_value(system_symbol, (yyvsp[(2) - (2)].dval));
          if (0 == status)
          {
            yyerror("System value specified on previous line is not valid");
            exit(1);
          }
        }
        else
        {
          status = dynsymbols_name_update((yyvsp[(1) - (2)].str), (yyvsp[(2) - (2)].dval), dynsymbols);
          if (0 == status)
          {
            yyerror("Error creating new parameter");
            assert(0);
            exit(1);
          }
        }
      }
    break;

  case 61:
#line 392 "parser.y"
    {
        int status;

        if (0 == strncasecmp((yyvsp[(1) - (3)].str), SYSTEM_PREFIX, sizeof(SYSTEM_PREFIX) - 1))
        {
          char *system_symbol;

          system_symbol = (yyvsp[(1) - (3)].str);

          // Skip over prefix.
          system_symbol += sizeof(SYSTEM_PREFIX) - 1;

          status = ds_system_set_value(system_symbol, (yyvsp[(3) - (3)].dval));
          if (0 == status)
          {
            yyerror("System value specified on previous line is not valid");
            exit(1);
          }
        }
        else
        {
          status = dynsymbols_name_update((yyvsp[(1) - (3)].str), (yyvsp[(3) - (3)].dval), dynsymbols);
          if (0 == status)
          {
            yyerror("Error creating new parameter");
            assert(0);
            exit(1);
          }
        }
      }
    break;

  case 63:
#line 430 "parser.y"
    { ; }
    break;

  case 64:
#line 435 "parser.y"
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
        rounded_pop = round((yyvsp[(3) - (3)].dval));
#else
        rounded_pop = floor((yyvsp[(3) - (3)].dval) + 0.5);
#endif

#ifndef CREATE_NON_AGGREGATE_SPECIE
        status = dynspecies_specie_add(current_specie_cgraph, (yyvsp[(1) - (3)].str), rounded_pop,
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
            status = dynspecies_specie_add(current_specie_cgraph, (yyvsp[(1) - (3)].str),
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
    break;

  case 65:
#line 500 "parser.y"
    {
        /* Finished initializing specie (always a single molecule). */
        initializing_specie = FALSE;

        /* Single molecule, so no inter-molecule connection to add. */

        /* Just a molecule declaration, so do not add to specie list. */
      }
    break;

  case 66:
#line 513 "parser.y"
    { ; }
    break;

  case 67:
#line 518 "parser.y"
    {
        /* Finished building rule. */
        building_rule = FALSE;
      }
    break;

  case 68:
#line 526 "parser.y"
    {
        int status;

        status = finish_reaction_rules((yyvsp[(1) - (5)].str), FALSE, (yyvsp[(2) - (5)].ival), (yyvsp[(4) - (5)].ival), (yyvsp[(5) - (5)].dval), 0.0);
        if (FALSE == status)
        {
          yyerror_prev("Error setting up reaction rule");
          assert(0);
          exit(1);
        }

      }
    break;

  case 69:
#line 539 "parser.y"
    {
        int status;

        status = finish_reaction_rules((yyvsp[(1) - (7)].str), TRUE, (yyvsp[(2) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(5) - (7)].dval), (yyvsp[(7) - (7)].dval));
        if (FALSE == status)
        {
          yyerror_prev("Error setting up reaction rule");
          exit(1);
        }
      }
    break;

  case 70:
#line 553 "parser.y"
    {
        /* Finished with left hand side. */
        reaction_lhs = FALSE;

        molecule_star_lhs = molecule_star_count;
        molecule_star_count = 0;

        /* Return size of list. */
        (yyval.ival) = (yyvsp[(1) - (1)].ival);
      }
    break;

  case 71:
#line 567 "parser.y"
    {
        /* Finished with right hand side. */

        molecule_star_rhs = molecule_star_count;
        molecule_star_count = 0;

        /* Return size of list. */
        (yyval.ival) = (yyvsp[(1) - (1)].ival);
      }
    break;

  case 72:
#line 580 "parser.y"
    {
        /* Keep track of connections for each specie. */
        building_rule = TRUE;

        (yyval.ival) = (yyvsp[(1) - (1)].ival);
      }
    break;

  case 73:
#line 587 "parser.y"
    {
        /* Make rule 'or' to lhs rule. */
        dynrules_rule_add_or_rule((yyvsp[(1) - (2)].ival), (yyvsp[(2) - (2)].ival), world->dynrules);

        (yyval.ival) = (yyvsp[(1) - (2)].ival);
      }
    break;

  case 74:
#line 597 "parser.y"
    {
        /* Specie is getting created or deleted. */
        molecule_star_count++;

        (yyval.ival) = 0;
      }
    break;

  case 75:
#line 604 "parser.y"
    {
        /* Keep track of connections for each specie. */
        building_rule = TRUE;

        /* Return size of list. */
        (yyval.ival) = 1;
      }
    break;

  case 76:
#line 612 "parser.y"
    {
        /* Specie is getting created or deleted. */
        molecule_star_count++;

        /* Return size of list. */
        (yyval.ival) = (yyvsp[(1) - (3)].ival);
      }
    break;

  case 77:
#line 620 "parser.y"
    {
        /* Return size of list. */
        (yyval.ival) = (yyvsp[(1) - (3)].ival) + 1;
      }
    break;

  case 78:
#line 628 "parser.y"
    {
        int status;

        building_rule = FALSE;

        // Passing a denominator of one means count each molecule.
        status = finish_observable((yyvsp[(3) - (4)].str), (yyvsp[(4) - (4)].ival), 1.0);
        if (FALSE == status)
        {
          yyerror_prev("Finishing molecule observable failed.");
          assert(0);
          exit(1);
        }
      }
    break;

  case 79:
#line 643 "parser.y"
    {
        int status;

        building_rule = FALSE;

        // Passing a denominator of zero means just count species,
        // no matter how many matching molecules found in each species.
        status = finish_observable((yyvsp[(3) - (4)].str), (yyvsp[(4) - (4)].ival), 0.0);
        if (FALSE == status)
        {
          yyerror_prev("Finishing specie observable failed.");
          assert(0);
          exit(1);
        }
      }
    break;

  case 80:
#line 662 "parser.y"
    {
        // Return last molecule rule index.
        (yyval.ival) = (yyvsp[(2) - (3)].ival);
      }
    break;

  case 81:
#line 667 "parser.y"
    {
        // Return last molecule rule index.
        (yyval.ival) = (yyvsp[(1) - (1)].ival);
      }
    break;

  case 82:
#line 672 "parser.y"
    {
        /* Make rule 'and' to lhs rule. */
        dynrules_rule_add_and_rule((yyvsp[(1) - (3)].ival), (yyvsp[(3) - (3)].ival), world->dynrules);

        (yyval.ival) = (yyvsp[(1) - (3)].ival);
      }
    break;

  case 83:
#line 679 "parser.y"
    {
        /* Make rule 'or' to lhs rule. */
        dynrules_rule_add_or_rule((yyvsp[(1) - (3)].ival), (yyvsp[(3) - (3)].ival), world->dynrules);

        // Return last molecule rule index.
        (yyval.ival) = (yyvsp[(1) - (3)].ival);
      }
    break;

  case 84:
#line 687 "parser.y"
    {
        /* Make rule 'xor' to lhs rule. */
        dynrules_rule_add_xor_rule((yyvsp[(1) - (3)].ival), (yyvsp[(3) - (3)].ival), world->dynrules);

        // Return last molecule rule index.
        (yyval.ival) = (yyvsp[(1) - (3)].ival);
      }
    break;

  case 85:
#line 698 "parser.y"
    {
        /* Finished building current specie. */
        building_specie = FALSE;

        // Return last molecule index.
        (yyval.ival) = (yyvsp[(1) - (1)].ival);
      }
    break;

  case 86:
#line 709 "parser.y"
    {
        /* Finished building current specie. */
        building_specie = FALSE;

        // Return last molecule index.
        (yyval.ival) = (yyvsp[(1) - (1)].ival);
      }
    break;

  case 87:
#line 720 "parser.y"
    {
        (yyval.ival) = (yyvsp[(1) - (1)].ival);
      }
    break;

  case 88:
#line 737 "parser.y"
    {
        (yyval.ival) = (yyvsp[(1) - (1)].ival);
      }
    break;

  case 89:
#line 741 "parser.y"
    {
        yyerror_prev("Modified molecule aggregate not yet supported");
        assert(0);
        exit(1);

        (yyval.ival) = (yyvsp[(2) - (4)].ival);
      }
    break;

  case 90:
#line 752 "parser.y"
    {
        /* Return lhs molecule index. */
        (yyval.ival) = (yyvsp[(1) - (1)].ival);
      }
    break;

  case 91:
#line 757 "parser.y"
    {
        // Don't add_and_rule because it breaks specie declaration for species with multiple molecules.
        //dynrules_rule_add_and_rule($1, $3, world->dynrules);
        /* Return lhs molecule index. */
        (yyval.ival) = (yyvsp[(1) - (3)].ival);
      }
    break;

  case 92:
#line 767 "parser.y"
    {
        /* Return lhs molecule index. */
        (yyval.ival) = (yyvsp[(1) - (1)].ival);
      }
    break;

  case 93:
#line 772 "parser.y"
    {
        dynrules_rule_add_and_rule((yyvsp[(1) - (3)].ival), (yyvsp[(3) - (3)].ival), world->dynrules);

        /* Return lhs molecule index. */
        (yyval.ival) = (yyvsp[(1) - (3)].ival);
      }
    break;

  case 94:
#line 782 "parser.y"
    {
        /* Finished with current molecule. */
        current_molecule = NULL;
        current_molecule_index = -1;
        current_molecule_rule_index = -1;
        current_object_index = -1;
        current_object_rule_index = -1;
        new_molecule = FALSE;

        /* Return molecule index. */
        (yyval.ival) = (yyvsp[(1) - (2)].ival);
      }
    break;

  case 95:
#line 798 "parser.y"
    {
        free(molecule_component_multipliers_table);
        molecule_component_multipliers_table = NULL;

        /* Return molecule index. */
        (yyval.ival) = (yyvsp[(1) - (2)].ival);
      }
    break;

  case 96:
#line 809 "parser.y"
    {
        free(molecule_component_multipliers_table);
        molecule_component_multipliers_table = NULL;

        /* Return molecule index. */
        (yyval.ival) = (yyvsp[(1) - (2)].ival);
      }
    break;

  case 97:
#line 820 "parser.y"
    {
        int molecule_color;
        boolean status;

        /* Create molecule to verify structure. */
        molecule_color = setup_molecule_name((yyvsp[(1) - (1)].str), FALSE);

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
        current_molecule_name = (yyvsp[(1) - (1)].str);

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
        (yyval.ival) = current_molecule_rule_index;
      }
    break;

  case 98:
#line 895 "parser.y"
    {
        int color;

        /* Create molecule to define structure or */
        /* verify structure if molecule already exists. */
        color = setup_molecule_name((yyvsp[(1) - (1)].str), TRUE);

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
        (yyval.ival) = current_molecule_index;
      }
    break;

  case 99:
#line 927 "parser.y"
    {
        /*
         * Now that components are finished, molecule modifiers will be
         * read if present.
         */
        current_object_index = current_molecule_index;
        current_object_name = current_molecule_name;
        current_object_rule_index = current_molecule_rule_index;
      }
    break;

  case 101:
#line 941 "parser.y"
    { ; }
    break;

  case 102:
#line 943 "parser.y"
    { ; }
    break;

  case 103:
#line 948 "parser.y"
    {
        boolean status;
        dynints_t *state_dynints;

        if (TRUE == initializing_specie)
        {
          state_dynints = current_specie_cgraph->
                          state_dynints_array.a[current_object_index];

          if (TRUE == new_molecule)
          {
            status = molecule_component_add((yyvsp[(1) - (2)].ival), (yyvsp[(2) - (2)].ival), state_dynints,
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
              molecule_component_verify((yyvsp[(1) - (2)].ival), (yyvsp[(2) - (2)].ival), state_dynints,
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
            molecule_component_verify((yyvsp[(1) - (2)].ival), (yyvsp[(2) - (2)].ival),
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
    break;

  case 104:
#line 1004 "parser.y"
    {
        color_t *component_color;

        /* Lookup component color. */
        component_color = real_dyncolors_color_create((yyvsp[(1) - (1)].str), COLOR_TYPE_COMPONENT,
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
          current_object_name = (yyvsp[(1) - (1)].str);

          /* Associate component with molecule. */
          dynrules_rule_add_adjacent_rule (current_molecule_rule_index,
                                           current_object_rule_index,
                                           world->dynrules);
        }

        (yyval.ival) = component_color->value;
      }
    break;

  case 105:
#line 1071 "parser.y"
    {
        (yyval.ival) = 1;
      }
    break;

  case 106:
#line 1075 "parser.y"
    {
        (yyval.ival) = (yyvsp[(1) - (1)].ival);
      }
    break;

  case 107:
#line 1082 "parser.y"
    {
        if (-1 == (yyvsp[(1) - (1)].ival))
        {
          (yyval.ival) = 1;
        }
        else
        {
          (yyval.ival) = (yyvsp[(1) - (1)].ival);
        }
      }
    break;

  case 108:
#line 1093 "parser.y"
    {
        if (-1 == (yyvsp[(1) - (2)].ival))
        {
          (yyval.ival) = (yyvsp[(2) - (2)].ival);
        }
        else if (-1 == (yyvsp[(2) - (2)].ival))
        {
          (yyval.ival) = (yyvsp[(1) - (2)].ival);
        }
        else
        {
          /* Two valid component multiples, so multiply them together. */
          (yyval.ival) = (yyvsp[(1) - (2)].ival) * (yyvsp[(2) - (2)].ival);
        }
      }
    break;

  case 109:
#line 1112 "parser.y"
    {
        (yyval.ival) = 1;
      }
    break;

  case 110:
#line 1116 "parser.y"
    {
        (yyval.ival) = (yyvsp[(1) - (1)].ival);
      }
    break;

  case 111:
#line 1123 "parser.y"
    {
        if (-1 == (yyvsp[(1) - (1)].ival))
        {
          (yyval.ival) = 1;
        }
        else
        {
          (yyval.ival) = (yyvsp[(1) - (1)].ival);
        }
      }
    break;

  case 112:
#line 1134 "parser.y"
    {
        if (-1 == (yyvsp[(1) - (2)].ival))
        {
          (yyval.ival) = (yyvsp[(2) - (2)].ival);
        }
        else if (-1 == (yyvsp[(2) - (2)].ival))
        {
          (yyval.ival) = (yyvsp[(1) - (2)].ival);
        }
        else
        {
          /* Two valid component multiples, so multiply them together. */
          (yyval.ival) = (yyvsp[(1) - (2)].ival) * (yyvsp[(2) - (2)].ival);
        }
      }
    break;

  case 113:
#line 1153 "parser.y"
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

        (yyval.ival) = -1;
      }
    break;

  case 114:
#line 1175 "parser.y"
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

        (yyval.ival) = -1;
      }
    break;

  case 115:
#line 1197 "parser.y"
    {
        int status;

        /* Add link name to component. */
        if (initializing_specie)
        {
          int status;

          /* Add current molecule component to link table */
          status = add_connection(specie_rule_indexes[current_specie_index],
                                  (yyvsp[(1) - (1)].str), current_object_index);
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
                                  (yyvsp[(1) - (1)].str), current_object_rule_index);
          if (FALSE == status)
          {
            yyerror_prev("Error adding link modifier");
            assert(0);
            exit(1);
          }
        }

        (yyval.ival) = -1;
      }
    break;

  case 116:
#line 1234 "parser.y"
    {
        (yyval.ival) = (yyvsp[(2) - (2)].ival);
      }
    break;

  case 117:
#line 1238 "parser.y"
    {
        (yyval.ival) = (yyvsp[(1) - (1)].ival);
      }
    break;

  case 118:
#line 1245 "parser.y"
    {
        int status;
        color_t *state_color;

        state_color = real_dyncolors_color_create((yyvsp[(2) - (2)].str), COLOR_TYPE_STATE,
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
        (yyval.ival) = -1;
      }
    break;

  case 119:
#line 1326 "parser.y"
    {
        yyerror_prev("Partitions not implemented yet");
        assert(0);
        exit(1);

        (yyval.ival) = -1;
      }
    break;

  case 120:
#line 1334 "parser.y"
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
          object_rule->count = (yyvsp[(2) - (2)].dval);
        }

        /* Pass multiple up the parse tree. */
        (yyval.ival) = (yyvsp[(2) - (2)].dval);
      }
    break;

  case 121:
#line 1358 "parser.y"
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
          object_rule->count = (yyvsp[(2) - (2)].dval) - 1;
        }

        (yyval.ival) = -1;
      }
    break;

  case 122:
#line 1380 "parser.y"
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
          object_rule->count = (yyvsp[(2) - (2)].dval) + 1;
        }

        (yyval.ival) = -1;
      }
    break;

  case 123:
#line 1402 "parser.y"
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
          object_rule->count = (yyvsp[(2) - (2)].dval);
        }

        (yyval.ival) = -1;
      }
    break;

  case 124:
#line 1424 "parser.y"
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
          object_rule->count = (yyvsp[(2) - (2)].dval);
        }

        (yyval.ival) = -1;
      }
    break;

  case 125:
#line 1449 "parser.y"
    { strcpy((yyval.str), (yyvsp[(1) - (1)].str)); }
    break;

  case 126:
#line 1451 "parser.y"
    { strcpy((yyval.str), (yyvsp[(1) - (1)].str)); }
    break;

  case 127:
#line 1456 "parser.y"
    { strcpy((yyval.str), (yyvsp[(1) - (1)].str)); }
    break;

  case 128:
#line 1458 "parser.y"
    { strcpy((yyval.str), (yyvsp[(1) - (1)].str)); }
    break;

  case 129:
#line 1460 "parser.y"
    { strcpy((yyval.str), "*"); }
    break;

  case 130:
#line 1465 "parser.y"
    {
        *(yyval.str) = '\0';
      }
    break;

  case 131:
#line 1469 "parser.y"
    {
        strcpy((yyval.str), (yyvsp[(1) - (1)].str));
      }
    break;

  case 132:
#line 1473 "parser.y"
    {
        strcpy((yyval.str), (yyvsp[(1) - (1)].str));
      }
    break;

  case 133:
#line 1477 "parser.y"
    {
        strcpy((yyval.str), (yyvsp[(1) - (1)].str));
      }
    break;

  case 134:
#line 1481 "parser.y"
    {
        strcpy((yyval.str), (yyvsp[(1) - (2)].str));
      }
    break;

  case 135:
#line 1488 "parser.y"
    {
        /* Finalize current hash. */
        ;
      }
    break;

  case 137:
#line 1497 "parser.y"
    {
        /* Create new hash structure and add current item to hash. */
        ;
      }
    break;

  case 138:
#line 1502 "parser.y"
    {
        /* Add current item to existing hash. */
        ;
      }
    break;

  case 139:
#line 1510 "parser.y"
    {
        /* Store current item. */
        yyerror_prev("Do not support hashes of hashes");
        assert(0);
        exit(1);
      }
    break;

  case 140:
#line 1517 "parser.y"
    {
        /* Store current item. */
        yyerror_prev("Do not support hashes of arrays");
        assert(0);
        exit(1);
      }
    break;

  case 141:
#line 1524 "parser.y"
    {
        /* Store current number. */
        bng_command_set_double((yyvsp[(1) - (3)].str), (yyvsp[(3) - (3)].dval));
      }
    break;

  case 142:
#line 1529 "parser.y"
    {
        /* Store current item. */
        bng_command_set_string((yyvsp[(1) - (3)].str), (yyvsp[(3) - (3)].str));
      }
    break;

  case 143:
#line 1537 "parser.y"
    {
        /* Store current item. */
        yyerror_prev("Do not support arrays of hashes");
        assert(0);
        exit(1);
      }
    break;

  case 144:
#line 1544 "parser.y"
    {
        /* Store current item. */
        yyerror_prev("Do not support arrays of arrays");
        assert(0);
        exit(1);
      }
    break;

  case 145:
#line 1551 "parser.y"
    {
        /* Store current number. */
        yyerror_prev("Do not support arrays of expressions");
        assert(0);
        exit(1);
      }
    break;

  case 146:
#line 1558 "parser.y"
    {
        /* Store current item. */
        yyerror_prev("Do not support arrays of strings");
        assert(0);
        exit(1);
      }
    break;

  case 147:
#line 1568 "parser.y"
    {
        /* Finalize current array. */
        ;
      }
    break;

  case 149:
#line 1577 "parser.y"
    { 
        ; 
      }
    break;

  case 150:
#line 1581 "parser.y"
    { 
        ; 
      }
    break;

  case 151:
#line 1587 "parser.y"
    { (yyval.dval) = (yyvsp[(1) - (1)].dval); }
    break;

  case 152:
#line 1588 "parser.y"
    { (yyval.dval) = (yyvsp[(1) - (3)].dval) + (yyvsp[(3) - (3)].dval); }
    break;

  case 153:
#line 1589 "parser.y"
    { (yyval.dval) = (yyvsp[(1) - (3)].dval) - (yyvsp[(3) - (3)].dval); }
    break;

  case 154:
#line 1590 "parser.y"
    { (yyval.dval) = (yyvsp[(1) - (3)].dval) * (yyvsp[(3) - (3)].dval); }
    break;

  case 155:
#line 1592 "parser.y"
    {
        if ((yyvsp[(3) - (3)].dval) == 0.0)
        {
          yyerror_prev("divide by zero");
          assert(0);
          exit(1);
        }
        else
        {
          (yyval.dval) = (yyvsp[(1) - (3)].dval) / (yyvsp[(3) - (3)].dval);
        }
      }
    break;

  case 156:
#line 1604 "parser.y"
    { (yyval.dval) = -(yyvsp[(2) - (2)].dval); }
    break;

  case 157:
#line 1605 "parser.y"
    { (yyval.dval) = pow((yyvsp[(1) - (3)].dval), (yyvsp[(3) - (3)].dval)); }
    break;

  case 158:
#line 1606 "parser.y"
    { (yyval.dval) = (yyvsp[(2) - (3)].dval); }
    break;

  case 159:
#line 1611 "parser.y"
    {
        char *endptr;

        (yyval.dval) = strtod((yyvsp[(1) - (1)].str), &endptr);
        if (endptr == (yyvsp[(1) - (1)].str))
        {
          yyerror_prev("Error occured when trying to convert string to number");
          assert(0);
          exit(1);
        }
      }
    break;

  case 160:
#line 1623 "parser.y"
    {
        size_t index;
        double value;

        if (0 == strncasecmp((yyvsp[(1) - (1)].str), SYSTEM_PREFIX, sizeof(SYSTEM_PREFIX) - 1))
        {
          char *system_symbol;

          system_symbol = (yyvsp[(1) - (1)].str);

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
          index = dynsymbols_get_value((yyvsp[(1) - (1)].str), &value, dynsymbols);
          if (0 == index)
          {
            yyerror("symbol not found");
            assert(0);
            exit(1);
          }
        }

        (yyval.dval) = value;
      }
    break;


/* Line 1267 of yacc.c.  */
#line 3471 "parser.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 1660 "parser.y"



