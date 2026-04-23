
/* The "incl" state is used for picking up the name of an include file. */
%x incl

%{
#include "parser.h"
#include "parser_lib.h"
#include "output.h"
#include <string.h>
#include <math.h>

#define YY_NO_UNPUT
#define MAX_INCLUDE_DEPTH 10

#define MYRETURN tokenpos += yyleng; return 

void suppress_warnings7() {SUPPRESS_WARNINGS;}

void yyerror(char *s);
void yyerror_prev(char *s);

/* For handling include files. */
YY_BUFFER_STATE include_stack[MAX_INCLUDE_DEPTH];
char *include_name_stack[MAX_INCLUDE_DEPTH];
int include_lineno_stack[MAX_INCLUDE_DEPTH];
int include_stack_index = 0;

/* For error reporting. */
char * current_filename = NULL; /* Don't have first filename */
static char linebuf[1024];
static int lineno = 1; /* First line will not be available for reporting. */
static int tokenpos = 0;
%}

qstring \"[^\"\n]*[\"\n]
sqstring '[^'\n]*['\n]

%%

        /* For error reporting. */
\n.* {
    /* Save the next line. */
    /* If line too long to fit in linebuf, everything will still work, */
    /* except if there is an error, the error message might look strange. */
    strncpy(linebuf, yytext+1, sizeof(linebuf));
    linebuf[sizeof(linebuf) - 1] = '\0';
    lineno++;
    tokenpos = 0;

    /* Give back all but the \n to rescan. */
    yyless(1);
  }

        /* For handling include files. */
include { BEGIN(incl); }
<incl>[ \t]* /* eat the whitespace */
<incl>[^ \t\n]+ { /* got the include file name */
    if (include_stack_index >= MAX_INCLUDE_DEPTH)
    {
      yyerror("Includes nested too deeply");
      exit(1);
    }

    include_stack[include_stack_index++] = YY_CURRENT_BUFFER;
    include_name_stack[include_stack_index] = current_filename;
    include_lineno_stack[include_stack_index] = lineno;

    if (current_filename != NULL)
    {
      free(current_filename);
    }

    yyin = fopen(yytext, "r");

    if (!yyin)
    {
      yyerror("Unable to open include file");
      exit(1);
    }

    current_filename = strdup(yytext);
    lineno = 1; /* First line will not be available for reporting. */
    linebuf[0] = '\0';

    yy_switch_to_buffer(yy_create_buffer(yyin, YY_BUF_SIZE));

    BEGIN(INITIAL);
  }

        /* Labels */
[ \t]*[A-Za-z0-9_-]*: {
    if (yyleng + 1 > sizeof(yylval.str))
    {
      yyerror("Label too long.");
      abort();
    }
    memcpy(yylval.str, yytext, yyleng + 1);

    /* Remove trailing colon. */
    yylval.str[yyleng - 2] = '\0';

    MYRETURN LABEL;
  }

	/* Block delimeters */
begin[ \t]+parameters { MYRETURN BEGIN_PARAMETERS; }
end[ \t]+parameters { MYRETURN END_PARAMETERS; }
begin[ \t]+species { MYRETURN BEGIN_SPECIES; }
begin[ \t]+seed[ \t]+species { MYRETURN BEGIN_SPECIES; }
end[ \t]+species { MYRETURN END_SPECIES; }
end[ \t]+seed[ \t]+species { MYRETURN END_SPECIES; }
begin[ \t]+reaction[ _]rules { MYRETURN BEGIN_REACTION_RULES; }
end[ \t]+reaction[ _]rules { MYRETURN END_REACTION_RULES; }
begin[ \t]+observables { MYRETURN BEGIN_OBSERVABLES; }
end[ \t]+observables { MYRETURN END_OBSERVABLES; }
begin[ \t]+molecule[ _]types { MYRETURN BEGIN_MOLECULE_TYPES;}
end[ \t]+molecule[ _]types { MYRETURN END_MOLECULE_TYPES;}
begin[ \t]+actions { MYRETURN BEGIN_ACTIONS; }
end[ \t]+actions { MYRETURN END_ACTIONS; }

	/* Reserved words */
molecules { MYRETURN MOLECULES; }
species { MYRETURN SPECIES; }

	/* Commands */
generate_network { MYRETURN GENERATE_NETWORK; }
simulate_ode { MYRETURN SIMULATE_ODE; }
simulate_ssa { MYRETURN SIMULATE_SSA; }
setup_ps { MYRETURN SETUP_PS; }
state_output_ps { MYRETURN STATE_OUTPUT_PS; }
simulate_ps { MYRETURN SIMULATE_PS; }
simulate_rs { MYRETURN SIMULATE_RS; }
simulate_rm { MYRETURN SIMULATE_RS; }
toSBMLfile { MYRETURN TO_SBML_FILE; }

	/* Logical operators */
and { MYRETURN LOGICAL_AND; }
or { MYRETURN LOGICAL_OR; }

	/* Punctuation */
"\?" { MYRETURN QUESTION_MARK; }
"->" { MYRETURN UNIDIRECTIONAL_RULE; }
"<->" { MYRETURN BIDIRECTIONAL_RULE; }
"=>" { MYRETURN ','; }  /* For perl hashes */
"<=" { MYRETURN CMP_LE; }
">=" { MYRETURN CMP_GE; }
"<" { MYRETURN CMP_LT; }
">" { MYRETURN CMP_GT; }
[-+*/(){};,.~!@=]  { MYRETURN yytext[0]; }

	/* Numbers */
(([0-9]+|([0-9]*\.[0-9]+))([eE][-+]?[0-9]+)?) {
    if (yyleng + 1 > sizeof(yylval.str))
    {
      yyerror("Number too long.");
      abort();
    }
    memcpy(yylval.str, yytext, yyleng + 1);

    MYRETURN NUMBER;
  }

	/* Names */
[A-Za-z][A-Za-z0-9_]* {
    if (yyleng + 1 > sizeof(yylval.str))
    {
      yyerror("Name too long.");
      abort();
    }
    memcpy(yylval.str, yytext, yyleng + 1);

    MYRETURN NAME;
  }

	/* Strings */
{qstring} {
    if (yyleng + 1 > sizeof(yylval.str))
    {
      yyerror("Quoted string too long.");
      abort();
    }
    memcpy(yylval.str, yytext + 1, yyleng + 1 - 1); /* Skip open quote. */

    if (yylval.str[yyleng - 2] != '"')
    {
      yyerror("Unterminated character string");
      abort();
    }
    else
    {
      /* Remove closing quote. */
      yylval.str[yyleng - 2] = '\0';
    }

    MYRETURN QSTRING;
  }

{sqstring} {
    if (yyleng + 1 > sizeof(yylval.str))
    {
      yyerror("Quoted string too long.");
      abort();
    }
    memcpy(yylval.str, yytext + 1, yyleng + 1 - 1); /* Skip open quote. */

    if (yylval.str[yyleng - 2] != '\'')
    {
      yyerror("Unterminated character string");
      abort();
    }
    else
    {
      /* Remove closing quote. */
      yylval.str[yyleng - 2] = '\0';
    }

    MYRETURN QSTRING;
  }

	/* Whitespace */
[ \t\r\n] {tokenpos += yyleng; /* Ignore whitespace. */}

\\ {tokenpos += yyleng; /* Ignore line continuation. */}

"#".*  {tokenpos += yyleng; /* Ignore comments. */}

<<EOF>> {
    if (--include_stack_index < 0)
    {
      yyterminate();
    }
    else
    {
      yy_delete_buffer(YY_CURRENT_BUFFER);
      yy_switch_to_buffer(include_stack[include_stack_index]);
      current_filename = include_name_stack[include_stack_index];
      lineno = include_lineno_stack[include_stack_index];
    }
  }

	/* Anything else */
.  { yyerror("Invalid character"); exit(1);}

%%

void yyerror(char *s)
{
  if ('\0' != linebuf[0])
  {
    if (NULL != current_filename)
    {
      error_printf("%s:", current_filename);
    }
    error_printf("%d: %s \n", lineno, s);
    error_printf("%s\n", linebuf);
    error_printf("%*s\n", 1 + tokenpos - yyleng, "^");
  }
  else
  {
    error_printf("line %d, column %d: %s \n",
               lineno, 1 + tokenpos - yyleng, s);
  }
}

void yyerror_prev(char *s)
{
  if (NULL != current_filename)
  {
    error_printf("%s:", current_filename);
  }
  error_printf("Somewhere on or above line %d: %s \n",
             lineno, s);
}

