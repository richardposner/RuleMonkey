## Makefile.am -- Process this ifile with automake to produce Makefile.in

AM_CFLAGS = @ac_CFLAGS@
#CFLAGS = -O0

noinst_DATA = cscope.out

cscope.out:
	if test -e '/usr/bin/cscope'; \
	then \
	 ls */*.[ch]* */*.[ly] > cscope.files; cscope -Rb | touch cscope.out; \
	fi

EXTRA_DIST   = COPYING.txt bootstrap models/A+A.bngl models/A+A_mixed_1.bngl \
	models/A+A_mixed_2.bngl models/A+B_rings.bngl \
	models/Tutorial_Example.bngl models/blbr.bngl models/egfr_net.bngl \
	models/fceri_ji.bngl models/lat.bngl models/mlnr.bngl \
	models/st.bngl models/st_multi_1.bngl models/st_multi_2.bngl \
	models/stiff.bngl models/testcase1.bngl models/testcase2.bngl \
	models/tlbr.bngl models/tlbr_rings.bngl models/toy-jim.bngl \
	models/pltr.bngl models/testcase2a.bngl models/testcase2b.bngl

SUBDIRS      = dSFMT nauty common particlestoc reactionstoc parser src

MAINTAINERCLEANFILES  = Makefile.in aclocal.m4 configure \
                        config-h.in nauty22/* dSFMT-src-1.3/*

