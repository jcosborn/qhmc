AM_CFLAGS = -I$(top_srcdir)/include -DLUA_COMPAT_ALL \
 $(QOPQDP_CFLAGS) \
 -DQOP_Colors="$(QOPQDP_COLORS)" \
 -DQDP_Colors="$(QOPQDP_COLORS)" \
 -DQLA_Colors="$(QOPQDP_COLORS)"

lib_LIBRARIES = ../libqhmc_qopqdp.a

HDRS = qhmc_qopqdp_common.h
SRC = \
 qhmc_qopqdp.c lattice.c subset.c \
 reader.c writer.c \
 lrstate.c \
 smear.c misc.c \
 asqtad.c asqtad-invert.c hisq.c squark.c \
 wilson.c wquark.c \
 dw.c dwquark.c
FSRC = \
 lrealF.c lcomplexF.c lcmatrixF.c \
 cscalarF.c gaugeF.c forceF.c gaugefixF.c
DSRC = \
 lrealD.c lcomplexD.c lcmatrixD.c \
 cscalarD.c gaugeD.c forceD.c gaugefixD.c

EXTRA_DIST = \
 lreal.c lcomplex.c lcmatrix.c \
 cscalar.c gauge.c force.c gaugefix.c

___libqhmc_qopqdp_a_SOURCES = \
 $(HDRS) $(SRC) $(FSRC) $(DSRC)
