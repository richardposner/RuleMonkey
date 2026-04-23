# Process this file with automake to produce Makefile.in

AM_CFLAGS = @ac_CFLAGS@

raw_source_dir = dSFMT
target_dir = $(top_builddir)/dSFMT-src-1.3

INCLUDES = -I$(target_dir)

EXTRA_DIST = dSFMT-src-1.3.tar.gz

noinst_DATA = $(target_dir)/.dSFMT_test

noinst_LTLIBRARIES  = libdsfmt.la
libdsfmt_la_SOURCES = \
	$(target_dir)/dSFMT.c \
	$(target_dir)/dSFMT.h \
	$(target_dir)/dSFMT-params.h \
	$(target_dir)/dSFMT-params11213.h \
	$(target_dir)/dSFMT-params1279.h \
	$(target_dir)/dSFMT-params132049.h \
	$(target_dir)/dSFMT-params19937.h \
	$(target_dir)/dSFMT-params216091.h \
	$(target_dir)/dSFMT-params2281.h \
	$(target_dir)/dSFMT-params4423.h \
	$(target_dir)/dSFMT-params44497.h \
	$(target_dir)/dSFMT-params607.h \
	$(target_dir)/dSFMT-params86243.h

# Don't do recursive libadd.
libdsfmt_la_LIBADD  =

$(target_dir)/.dSFMT_test: dSFMT-src-1.3.tar.gz
	cd $(top_builddir) && \
	gunzip -c $(raw_source_dir)/dSFMT-src-1.3.tar.gz | \
	tar -xf - && \
	cd - && \
	touch $(target_dir)/.dSFMT_test

$(libdsfmt_la_SOURCES): $(target_dir)/.dSFMT_test
