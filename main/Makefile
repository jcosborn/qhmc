include ../make.inc

DEFS =
ifneq ($(strip $(HAVE_QOPQDP)),)
  include ../qopqdp/make-qopqdp.inc
  DEFS += -DHAVE_QOPQDP
endif

CFLAGS = $(COPT) $(DEFS) $(LUAINC)
LDFLAGS = $(COPT)
LDLIBS = $(QOPQDPLIBS) $(LUALIBS)
DEPS = $(QOPQDPDEPS) $(LUADEPS)

PROGS = qhmc
OBJS = lua.o qhmc.o

all: $(PROGS)
$(PROGS): $(OBJS) $(DEPS)
$(OBJS): qhmc.h

clean:
	-rm -f $(OBJS) $(PROGS)

realclean: clean
