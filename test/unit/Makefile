include ../../make.inc
include ../../qopqdp/make-qopqdp.inc

CFLAGS = $(COPT) $(QOPQDPINC) $(LUAINC) -I../../qopqdp
LDFLAGS = $(COPT)
LDLIBS = $(LQOP) $(LQDP) $(LQLA) $(LQIO) $(LQMP) -lm

HDRS = ../../qopqdp/qhmc_qopqdp.h ../../qopqdp/qhmc_qopqdp_common.h
OBJS = ../../qopqdp/misc.o

all: sqrt

sqrt: $(OBJS)

clean:
	-rm -f sqrt $(OBJS)

realclean: clean
