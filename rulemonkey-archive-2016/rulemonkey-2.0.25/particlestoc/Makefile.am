# Process this file with automake to produce Makefile.in

MAINTAINERCLEANFILES = Makefile.in

AM_CFLAGS = @ac_CFLAGS@ @ac_CFLAGS_W@

INCLUDES = @ac_INCLUDES@

noinst_LTLIBRARIES    = libparticlestoc.la
libparticlestoc_la_SOURCES = \
	dynps_reactant_containers.c \
	dynps_reactant_containers.h \
	dynps_reaction_containers.c \
	dynps_reaction_containers.h \
	dynps_species.c \
	dynps_species.h \
	particlestoc.c \
	particlestoc.h \
	ps_particle.c \
	ps_particle.h \
	ps_reactant_container.c \
	ps_reactant_container.h \
	ps_reaction_container.c \
	ps_reaction_container.h \
	ps_reaction_definition.c \
	ps_reaction_definition.h \
	ps_specie.c \
	ps_specie.h

# Don't do recursive libadd.
libparticlestoc_la_LIBADD  =

