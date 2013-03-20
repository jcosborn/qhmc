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

SUBDIRS = $(DIRS) main

.PHONY: all subdirs clean realclean lua luaclean $(SUBDIRS)

all: lua subdirs
main: lua $(DIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean: $(SUBDIRS)

realclean: $(SUBDIRS) luaclean

lua:
	$(MAKE) -C $(LUADIR) PLAT=generic CC="$(CC)" CFLAGS="$(COPT) $(LUAOPTS)"

luaclean:
	$(MAKE) -C $(LUADIR) clean

VER = $(shell date +%Y%m%d%H%M)
TARNAME = qhmc-$(VER).tar.gz
tar:
	touch $(TARNAME)
	tar zcvf $(TARNAME) --transform 's|^./|qhmc-$(VER)/|' -X excludes.txt .
