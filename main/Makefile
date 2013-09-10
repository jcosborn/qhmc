include ../make.inc

DEFS =
ifneq ($(strip $(LFSDIR)),)
  include ../$(LFSDIR)/make-lfs.inc
  DEFS += -DHAVE_LFS
endif
ifneq ($(strip $(HAVE_QOPQDP)),)
  include ../qopqdp/make-qopqdp.inc
  DEFS += -DHAVE_QOPQDP
endif

CFLAGS = $(COPT) $(DEFS) $(LFSINC) $(LUAINC)
LDFLAGS = $(COPT)
LDLIBS = $(QOPQDPLIBS) $(LFSOBJS) $(LUALIBS) ../lbc/liblbc.a
DEPS = $(QOPQDPDEPS) $(LFSDEPS) $(LUADEPS)

PROGS = qhmc
OBJS = lua.o qhmc.o qhmc_common.o qhmc_complex.o

all: $(PROGS)
$(PROGS): $(OBJS) $(DEPS)
$(OBJS): qhmc.h qhmc_common.h
qhmc_complex.o: qhmc_complex.h

clean:
	-rm -f $(OBJS) $(PROGS)

realclean: clean
