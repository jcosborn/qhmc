include ../../make.inc
include ../../qopqdp/make-qopqdp.inc

LUAINC = -I../../$(LUADIR)/src -DLUA_COMPAT_ALL
CFLAGS = $(COPT) $(QOPQDPINC) $(LUAINC) -I../../qopqdp
LDFLAGS = $(COPT)
LDLIBS = $(LQOP) $(LQDP) $(LQLA) $(LQIO) $(LQMP) -lm

HDRS = ../../qopqdp/qhmc_qopqdp.h ../../qopqdp/qhmc_qopqdp_common.h
OBJS = ../../qopqdp/misc.o

all: exp

exp: $(OBJS)

clean:
	-rm -f exp $(OBJS)

realclean: clean
