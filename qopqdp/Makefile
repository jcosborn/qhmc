include ../make.inc
include make-qopqdp.inc

CFLAGS0 = $(COPT) $(QOPQDPINC) $(LUAINC) -I../main
FFLAGS = -DQOP_Precision=\'F\' -DQDP_Precision=\'F\' -DQLA_Precision=\'F\'
DFLAGS = -DQOP_Precision=\'D\' -DQDP_Precision=\'D\' -DQLA_Precision=\'D\'
CFLAGS = $(CFLAGS0) $(DFLAGS)

HDRS = qhmc_qopqdp.h qhmc_qopqdp_common.h
OBJS0 = \
 qhmc_qopqdp.o lattice.o \
 reader.o writer.o \
 smear.o misc.o \
 asqtad.o asqtad-invert.o hisq.o squark.o \
 wilson.o wquark.o \
 dw.o dwquark.o
FDOBJS = \
 cscalar.o \
 gauge.o force.o gaugefix.o
FOBJS = $(patsubst %.o,%F.o,$(FDOBJS))
DOBJS = $(patsubst %.o,%D.o,$(FDOBJS))
OBJS = $(OBJS0) $(FOBJS) $(DOBJS)

all: $(QHMCQOPQDP)

$(QHMCQOPQDP): $(OBJS)
	ar cru $@ $^

$(OBJS): $(HDRS)

$(FOBJS): %F.o:%.c
	$(CC) $(CFLAGS0) $(FFLAGS) -c $< -o $@

$(DOBJS): %D.o:%.c
	$(CC) $(CFLAGS0) $(DFLAGS) -c $< -o $@

clean:
	-rm -f $(QHMCQOPQDP) $(OBJS)

realclean: clean
