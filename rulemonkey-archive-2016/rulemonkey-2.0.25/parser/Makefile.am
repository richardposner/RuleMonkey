# Process this file with automake to produce Makefile.in

AM_CFLAGS = @ac_CFLAGS@ @ac_CFLAGS_W@

MAINTAINERCLEANFILES = Makefile.in

INCLUDES	= @ac_INCLUDES@

noinst_LTLIBRARIES = libparser.la

BUILT_SOURCES = parser.h
AM_YFLAGS = -d $(ac_YFLAGS)
YYDEBUG = @ac_YYDEBUG@
AM_LFLAGS = -i $(ac_LFLAGS)
libparser_la_SOURCES = parser.y lexer.l parser_lib.c parser_lib.h

libparser_la_LIBADD  = \
	$(top_builddir)/nauty/libnauty.la \
	$(top_builddir)/dSFMT/libdsfmt.la \
	$(top_builddir)/common/libcommon.la \
	$(top_builddir)/particlestoc/libparticlestoc.la \
	$(top_builddir)/reactionstoc/libreactionstoc.la \
	@LEXLIB@

