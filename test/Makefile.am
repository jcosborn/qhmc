#EXTRA_DIST =

#TESTS_ENVIRONMENT = @BACKEND_TEST_CMD@
#TESTS_ENVIRONMENT = $(top_srcdir)/test/runtest $(top_builddir)/bin/qhmc

TEST_EXTENSIONS = .lua
LUA_LOG_COMPILER = $(srcdir)/runtest $(top_builddir)/bin/qhmc

TESTS = \
 subset.lua \
 norm2.lua \
 point.lua \
 momentum.lua \
 fieldIO.lua \
 milc-su3_rhmc_hisq.lua \
 wilsonHMC.lua \
 wilson2fHMC.lua \
 hyp8fHMC.lua

EXTRA_DIST = \
 runtest runtest.lua \
 subset.lua subset.ref \
 norm2.lua norm2.ref \
 momentum.lua momentum.ref \
 point.lua point.ref \
 fieldIO.lua fieldIO.ref \
 milc-su3_rhmc_hisq.lua milc-su3_rhmc_hisq.in milc-su3_rhmc_hisq.rat milc-su3_rhmc_hisq.ref \
 wilsonHMC.lua wilsonHMC.ref \
 wilson2fHMC.lua wilson2fHMC.ref \
 hyp8fHMC.lua hyp8fHMC.ref

distclean-local:
	-rm -f *.out
