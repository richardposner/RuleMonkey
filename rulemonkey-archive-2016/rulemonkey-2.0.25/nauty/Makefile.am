# Process this file with automake to produce Makefile.in

AM_CFLAGS = @ac_CFLAGS@

raw_source_dir = nauty
target_dir = $(top_builddir)/nauty22

INCLUDES = -I$(target_dir)

EXTRA_DIST = nauty22.tar.gz

noinst_LTLIBRARIES  = libnauty.la
libnauty_la_SOURCES = \
	$(target_dir)/nautinv.c \
	$(target_dir)/naututil.c \
	$(target_dir)/nauty.c \
	$(target_dir)/rng.c \
	$(target_dir)/naugraph.c \
	$(target_dir)/nautil.c \
	$(target_dir)/nautinv.h \
	$(target_dir)/naututil.h \
	$(target_dir)/nauty.h \
	$(target_dir)/rng.h

noinst_DATA = $(target_dir)/config.status $(libnauty_la_SOURCES)

# Don't do recursive libadd.
libnauty_la_LIBADD  =

$(target_dir)/configure: nauty22.tar.gz
	cd $(top_builddir) && \
	gunzip -c $(raw_source_dir)/nauty22.tar.gz | \
	tar -xf - && \
	cd - && \
	touch $(target_dir)/configure

$(target_dir)/config.status: $(target_dir)/configure
	cd $(target_dir)/ && ./configure && cd -

$(libnauty_la_SOURCES): $(target_dir)/config.status

