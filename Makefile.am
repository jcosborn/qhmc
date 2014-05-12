SUBDIRS = . include lib bin test .

DISTCHECK_CONFIGURE_FLAGS = @DISTCHECK_CONFIGURE_FLAGS@

EXTRA_DIST = examples hmc

#all-local: timestamp $(abs_top_builddir)/hmc
all-local: timestamp

timestamp:
	@-echo "timestamp: `date`"

#$(abs_top_builddir)/hmc:
#	$(LN_S) $(top_srcdir)/hmc $@

uninstall-local:
	if [ -e $(DESTDIR)$(prefix)/hmc ]; then \
	  chmod -R 755 $(DESTDIR)$(prefix)/hmc; \
	  rm -rf $(DESTDIR)$(prefix)/hmc; \
	fi

install-data-local: uninstall-local
	if [ -e $(DESTDIR)$(prefix) ]; then \
	  cp -a $(top_srcdir)/hmc $(DESTDIR)$(prefix); \
	fi
