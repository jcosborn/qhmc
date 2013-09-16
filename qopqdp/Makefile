include ../make.inc
include make-qopqdp.inc

CFLAGS = $(COPT) $(QOPQDPINC) $(LUAINC) -I../main

HDRS = qhmc_qopqdp.h qhmc_qopqdp_common.h
OBJS = qhmc_qopqdp.o misc.o reader.o writer.o \
 cscalar.o \
 gauge.o gaugefix.o force.o smear.o \
 asqtad.o asqtad-invert.o hisq.o squark.o \
 wilson.o wquark.o \
 dw.o dwquark.o

all: $(QHMCQOPQDP)

$(QHMCQOPQDP): $(OBJS)
	ar cru $@ $^

$(OBJS): $(HDRS)

clean:
	-rm -f $(QHMCQOPQDP) $(OBJS)

realclean: clean
