#EXTRA_DIST =

#TESTS_ENVIRONMENT = @BACKEND_TEST_CMD@
#TESTS_ENVIRONMENT = $(top_srcdir)/test/runtest $(top_builddir)/bin/qhmc

TEST_EXTENSIONS = .lua
LUA_LOG_COMPILER = $(srcdir)/runtest $(top_builddir)/bin/qhmc

TESTS = \
 fieldIO.lua \
 gauge.lua \
 momentum.lua \
 norm2.lua \
 point.lua \
 random.lua \
 subset.lua \
 milc-su3_rhmc_hisq.lua \
 wilsonHMC.lua \
 wilson2fHMC.lua \
 hyp8fHMC.lua

EXTRA_DIST = \
 runtest runtest.lua \
 fieldIO.lua fieldIO.raw fieldIO.ref \
 gauge.lua gauge.raw gauge.ref \
 momentum.lua momentum.raw momentum.ref \
 norm2.lua norm2.raw norm2.ref \
 point.lua point.raw point.ref \
 random.lua random.raw random.ref \
 subset.lua subset.raw subset.ref \
 milc-su3_rhmc_hisq.lua milc-su3_rhmc_hisq.raw milc-su3_rhmc_hisq.ref \
  milc-su3_rhmc_hisq.in milc-su3_rhmc_hisq.rat \
 wilsonHMC.lua wilsonHMC.raw wilsonHMC.ref \
 wilson2fHMC.lua wilson2fHMC.raw wilson2fHMC.ref \
 hyp8fHMC.lua hyp8fHMC.raw hyp8fHMC.ref

distclean-local:
	-rm -f *.out *.tmp
