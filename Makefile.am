SUBDIRS = . include lib bin .

DISTCHECK_CONFIGURE_FLAGS = @DISTCHECK_CONFIGURE_FLAGS@

EXTRA_DIST = examples hmc

all-local: timestamp $(abs_top_builddir)/hmc

timestamp:
	@-echo "timestamp: `date`"

$(abs_top_builddir)/hmc:
	$(LN_S) $(top_srcdir)/hmc $@

install-data-local:
	\cp -a $(top_srcdir)/hmc $(prefix)
