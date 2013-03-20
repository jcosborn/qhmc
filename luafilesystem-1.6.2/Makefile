include ../make.inc
include make-lfs.inc
ifneq ($(strip $(HAVE_QOPQDP)),)
  include ../qopqdp/make-qopqdp.inc
endif

CFLAGS = $(COPT) $(LUAINC)

HDRS = src/lfs.h
OBJS = $(LFSOBJS)

all: $(OBJS)

$(OBJS): $(HDRS)

clean:
	-rm -f $(OBJS)

realclean: clean
