# Process this file with automake to produce Makefile.in

#AM_CFLAGS = @ac_CFLAGS@ @ac_CFLAGS_W@
AM_CFLAGS = @ac_CFLAGS@ -DPACKAGE_STRING="\"@PACKAGE_STRING@\""
AM_LDFLAGS = @ac_LDFLAGS@

MAINTAINERCLEANFILES = Makefile.in

INCLUDES	= @ac_INCLUDES@

bin_PROGRAMS = rm_bngl
rm_bngl_SOURCES = rm_bngl.c

rm_bngl_LDADD  = $(top_builddir)/parser/libparser.la


