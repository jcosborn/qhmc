include ../make.inc
include make-qopqdp.inc

CFLAGS = $(COPT) $(QOPQDPINC) $(LUAINC)

HDRS = qhmc_qopqdp.h qhmc_qopqdp_common.h
OBJS = qhmc_qopqdp.o misc.o gauge.o force.o smear.o \
 hisq.o hisq-invert.o squark.o \
 wilson.o wquark.o \
 dw.o dwquark.o

all: $(QHMCQOPQDP)

$(QHMCQOPQDP): $(OBJS)
	ar cru $@ $^

$(OBJS): $(HDRS)

clean:
	-rm -f $(QHMCQOPQDP) $(OBJS)

realclean: clean
