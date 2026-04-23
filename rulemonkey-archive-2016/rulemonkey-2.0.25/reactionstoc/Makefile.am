# Process this file with automake to produce Makefile.in

MAINTAINERCLEANFILES = Makefile.in

AM_CFLAGS = @ac_CFLAGS@ @ac_CFLAGS_W@

INCLUDES = @ac_INCLUDES@

noinst_LTLIBRARIES    = libreactionstoc.la
libreactionstoc_la_SOURCES = \
	reactionstoc.c \
	reactionstoc.h \
	rs_reaction_definition.c \
	rs_reaction_definition.h \
	rs_specie.c \
	rs_specie.h \
	rs_specie_container.c \
	rs_specie_container.h

# Don't do recursive libadd.
libreactionstoc_la_LIBADD  =

