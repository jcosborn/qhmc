SUBDIRS = . include lib bin test .

DISTCHECK_CONFIGURE_FLAGS = @DISTCHECK_CONFIGURE_FLAGS@

EXTRA_DIST = examples lua

#all-local: timestamp $(abs_top_builddir)/lua
all-local: timestamp

timestamp:
	@-echo "timestamp: `date`"

#$(abs_top_builddir)/lua:
#	$(LN_S) $(top_srcdir)/lua $@

uninstall-local:
	if [ -e $(DESTDIR)$(prefix)/lua ]; then \
	  chmod -R 755 $(DESTDIR)$(prefix)/lua; \
	  rm -rf $(DESTDIR)$(prefix)/lua; \
	  rm -rf $(DESTDIR)$(prefix)/hmc; \
	fi

install-data-local: uninstall-local
	if [ -e $(DESTDIR)$(prefix) ]; then \
	  cp -a $(top_srcdir)/lua $(DESTDIR)$(prefix); \
	  ln -s ../lua $(DESTDIR)$(prefix)/hmc; \
	fi
