#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

AC_INIT([rulemonkey], [2.0.25], [rulemonkey@tgen.org])
AM_CONFIG_HEADER([config/config.h])
AC_CONFIG_SRCDIR([common/world.c])
AC_CONFIG_AUX_DIR([config])

AM_INIT_AUTOMAKE([1.9 foreign])

cflags_optimize="-O3"
cflags_warnings="-Wall -Wextra -Wwrite-strings -Wconversion -Wfloat-equal -Wtradional -Wshadow -Wpointer-arith -Wshorten-64-to-32 -Wsign-compare -Wmissing-prototypes -Wmissing-field-initializers -Wunreachable-code -Wconversion"

ac_INCLUDES="\
 -I../parser \
 -I../common \
 -I../particlestoc \
 -I../reactionstoc \
 -I../nauty22 \
 -I../dSFMT-src-1.3"

# Consider -DNDEBUG when no double checking.

cflags_g4="-malign-natural -fast -mcpu=G4 -mtune=G4 -mpowerpc-gpopt -fgcse -mdynamic-no-pic -falign-loops=16"
cflags_g5="-malign-natural -fast -mcpu=G5 -mtune=G5 -mpowerpc-gpopt -fgcse -mdynamic-no-pic -mpowerpc64 -falign-loops=16"
cflags_prescott="-fast -march=prescott -mfpmath=sse,387 -mieee-fp -ftree-vectorize -fivopts -ftree-loop-linear"
cflags_intel64="-fast -m64 -mfpmath=sse,387 -mieee-fp -ftree-vectorize -fivopts -ftree-loop-linear"
cglags_64="-m64"

cflags_output_flush="-DOUTPUT_FLUSH"

# Always define DSFMT_MEXP to avoid compiler warnings.
# Always enable mersene twister PRNG.
cflags_mt_common="-DDSFMT_MEXP=19937"

cflags_mt_gcc="-DRANDOM_MT"
cflags_mt_gcc_sse="-DHAVE_SSE2 -msse2"
cflags_mt_gcc4_extra="--param max-inline-insns-single=1800 --param inline-unit-growth=500 --param large-function-growth=900"



# Check for parameters.
AC_ARG_ENABLE(debug,AS_HELP_STRING([--enable-debug],[enable additional debugging code]),,enable_debug=no)
AC_ARG_ENABLE(opt_debug,AS_HELP_STRING([--enable-opt-debug],[enable optimized debugging to find uninitialized variable usage]),,enable_opt_debug=no)
AC_ARG_ENABLE(yacc_debug,AS_HELP_STRING([--enable-yacc-debug],[enable yacc debugging to print extra yacc information]),,enable_yacc_debug=no)
AC_ARG_ENABLE(lex_debug,AS_HELP_STRING([--enable-lex-debug],[enable lex debugging to print extra lex information]),,enable_lex_debug=no)
AC_ARG_ENABLE(rule_create_verbose,AS_HELP_STRING([--enable-rule-create-verbose],[enable extra output related to creating rules]),,enable_rule_create_verbose=no)
AC_ARG_ENABLE(debug,AS_HELP_STRING([--enable-aggressive-check],[enable additional aggressive checking code]),,enable_aggressive_check=no)
AC_ARG_ENABLE(rule_check_verbose,AS_HELP_STRING(--enable-rule-check-verbose,enable extra output related to checking rules),,enable_rule_check_verbose=no)
AC_ARG_ENABLE(random_mt,AS_HELP_STRING(--disable-random-mt, Do not use Mersenne Twister random number generator),,enable_random_mt=yes)
AC_ARG_ENABLE(random_mt_sse,AS_HELP_STRING(--enable-random-mt-sse, Use SSE accellerated Mersenne Twister random number generator),,enable_random_mt_sse=no)
AC_ARG_ENABLE(cflags_w,AS_HELP_STRING(--enable-cflags-w,enable including -Wall -Werror),,enable_cflags_w=no)
AC_ARG_ENABLE(cflags_w,AS_HELP_STRING(--enable-output-flush,flush streams every time something is output),,enable_output_flush=no)
AC_ARG_ENABLE(intel64,AS_HELP_STRING(--enable-intel64,enable 64 bit Intel support),,enable_intel64=no)
AC_ARG_ENABLE(64,AS_HELP_STRING(--enable-64,enable 64 bit Intel support),,enable_64=no)

# Add compiler switches as needed.
cflags_debug=""
cflags_verbose=""
if test "$enable_debug" = "yes"; then
	cflags_debug="$cflags_debug -g"
fi

if test "$enable_opt_debug" = "yes"; then
	cflags_debug="$cflags_debug $cflags_optimize -g"
fi

if test "$enable_aggressive_check" = "yes"; then
	cflags_verbose="$cflags_verbose -DAGGRESSIVE_CHECK"
fi

if test "$enable_rule_create_verbose" = "yes"; then
	cflags_verbose="$cflags_verbose -DRULE_CREATE_DEBUG"
fi

if test "$enable_rule_check_verbose" = "yes"; then
  	cflags_verbose="$cflags_verbose -DRULE_CHECK_DEBUG"
fi

orig_cflags=$CFLAGS

# Checks for programs.
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_YACC
AC_PROG_LEX

CFLAGS=$orig_cflags

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([limits.h malloc.h stddef.h stdlib.h string.h strings.h sys/time.h unistd.h stdio.h assert.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_HEADER_TIME

#AC_SUBST(AM_CFLAGS)

orig_cflags=$CFLAGS
cflags_arch=""
cflags_mt=""
CFLAGS="$cflags_g4"
AC_RUN_IFELSE([AC_LANG_PROGRAM()], cflags_arch="$CFLAGS")
CFLAGS="$cflags_g5"
AC_RUN_IFELSE([AC_LANG_PROGRAM()], cflags_arch="$CFLAGS")
CFLAGS="$cflags_prescott"
AC_RUN_IFELSE([AC_LANG_PROGRAM()], cflags_arch="$CFLAGS")
CFLAGS="$cflags_mt_common"
AC_RUN_IFELSE([AC_LANG_PROGRAM()], cflags_mt="$CFLAGS")
CFLAGS="$cflags_mt_common $cflags_mt_gcc"
AC_RUN_IFELSE([AC_LANG_PROGRAM()], cflags_mt="$CFLAGS")
CFLAGS="$cflags_mt_common $cflags_mt_gcc $cflags_mt_gcc4_extra"
AC_RUN_IFELSE([AC_LANG_PROGRAM()], cflags_mt="$CFLAGS")
CFLAGS=$orig_cflags

if test "$enable_intel64" = "yes"; then
  	cflags_arch="$cflags_intel64 $cflags_mt_sse $cflags_mt"
	PROJECT_STRING="$PROJECT_STRING 64 bit";
fi
if test "$enable_64" = "yes"; then
  	cflags_arch="$cflags_64"
	PROJECT_STRING="$PROJECT_STRING 64 bit";
fi

ac_LDFLAGS="-lm"
ac_CFLAGS="$CFLAGS $cflags_mt_common -Wall"
ac_CFLAGS_W=""
ac_YFLAGS=""
ac_LFLAGS=""
ac_YYDEBUG=0

if test "$enable_cflags_w" = "yes"; then
	ac_CFLAGS_W="-Werror"
fi

if test "$enable_output_flush" = "yes"; then
  	ac_CFLAGS="$ac_CFLAGS $cflags_output_flush"
fi

if test "$enable_random_mt" = "yes"; then
  	ac_CFLAGS="$ac_CFLAGS $cflags_mt"
fi

if test "$enable_random_mt_sse" = "yes"; then
  	ac_CFLAGS="$ac_CFLAGS $cflags_mt_sse $cflags_mt"
fi

if test "$enable_debug" = "yes"; then
	ac_CFLAGS="$ac_CFLAGS $cflags_verbose $cflags_debug"
else
	# The variable cflags_debug will be set if --enable_opt_debug was given.
	ac_CFLAGS="$ac_CFLAGS $cflags_verbose $cflags_debug $cflags_arch $cflags_optimize"
fi

if test "$enable_yacc_debug" = "yes"; then
	ac_YYDEBUG=1
	ac_YFLAGS="$ac_YFLAGS -t -v"
fi

if test "$enable_lex_debug" = "yes"; then
	ac_LFLAGS="-d"
fi

AC_SUBST(ac_YYDEBUG)
AC_SUBST(ac_INCLUDES)
AC_SUBST(ac_YFLAGS)
AC_SUBST(ac_LFLAGS)
AC_SUBST(ac_LDFLAGS)
AC_SUBST(ac_CFLAGS)
AC_SUBST(ac_CFLAGS_W)
AC_SUBST(cflags_debug)
AC_SUBST(enable_debug)

# Checks for library functions.
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_TYPE_SIGNAL
AC_FUNC_STRTOD
AC_FUNC_ALLOCA
AC_CHECK_FUNCS([dup2 gettimeofday memset strcasecmp strdup round asprintf])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 parser/Makefile
                 common/Makefile
                 particlestoc/Makefile
                 reactionstoc/Makefile
                 dSFMT/Makefile
                 nauty/Makefile])

AC_OUTPUT

