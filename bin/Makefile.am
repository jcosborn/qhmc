AM_CFLAGS = -I$(top_srcdir)/include

QLIBS = ../lib/libqhmc_base.a
XLIBS = 
if HAVEQUDA
  QLIBS += ../lib/libqhmc_quda.a
  XLIBS += $(QUDA_LIBS)
endif
if HAVEQOPQDP
  QLIBS += ../lib/libqhmc_qopqdp.a
  XLIBS += $(QOPQDP_LIBS)
endif
if HAVELFS
  QLIBS += ../lib/liblfs.a
endif
QLIBS += ../lib/liblbc.a
QLIBS += ../lib/libqhmc_base.a
QLIBS += ../lib/liblua.a
LDADD = $(QLIBS) $(XLIBS)

bin_PROGRAMS = qhmc

qhmc_SOURCES = \
 qhmc_main.h lua.c qhmc.c
