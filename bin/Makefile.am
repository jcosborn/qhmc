AM_CFLAGS = -I$(top_srcdir)/include

LIBS += ../lib/libqhmc_base.a
if HAVEQUDA
  LIBS += ../lib/libqhmc_quda.a $(QUDA_LIBS)
endif
if HAVEQOPQDP
  LIBS += ../lib/libqhmc_qopqdp.a $(QOPQDP_LIBS)
endif
if HAVELFS
  LIBS += ../lib/liblfs.a
endif
LIBS += ../lib/liblbc.a
LIBS += ../lib/libqhmc_base.a
LIBS += ../lib/liblua.a

bin_PROGRAMS = qhmc

qhmc_SOURCES = \
 qhmc_main.h lua.c qhmc.c
