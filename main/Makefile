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
LDLIBS = $(QOPQDPLIBS) $(LFSOBJS) $(LUALIBS)
DEPS = $(QOPQDPDEPS) $(LFSDEPS) $(LUADEPS)

PROGS = qhmc
OBJS = lua.o qhmc.o

all: $(PROGS)
$(PROGS): $(OBJS) $(DEPS)
$(OBJS): qhmc.h

clean:
	-rm -f $(OBJS) $(PROGS)

realclean: clean
