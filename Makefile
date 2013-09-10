include make.inc

DIRS =
ifneq ($(strip $(LFSDIR)),)
  DIRS += $(LFSDIR)
  include $(LFSDIR)/make-lfs.inc
endif
ifneq ($(strip $(HAVE_QOPQDP)),)
  DIRS += qopqdp
  include qopqdp/make-qopqdp.inc
endif
SUBDIRS = $(DIRS) lbc main

.PHONY: all clean realclean lua $(SUBDIRS)

all: lua $(SUBDIRS)
main: lua $(DIRS)  # needed for parallel make to esure main is last
$(LFSDIR): lua
clean: $(SUBDIRS)
realclean: $(SUBDIRS) lua

ifeq ($(MAKECMDGOALS),realclean)
  lua: luaclean
else
  ifneq ($(MAKECMDGOALS),clean)
    lua: luabuild
  endif
endif

luabuild:
	$(MAKE) -C $(LUADIR) PLAT=generic CC="$(CC)" CFLAGS="$(COPT) $(LUAOPTS)"

luaclean:
	$(MAKE) -C $(LUADIR) clean

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

VER = $(shell date +%Y%m%d%H%M)
TARNAME = qhmc-$(VER).tar.gz
tar:
	touch $(TARNAME)
	tar zcvf $(TARNAME) --transform 's|^./|qhmc-$(VER)/|' -X excludes.txt .
