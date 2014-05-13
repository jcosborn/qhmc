#EXTRA_DIST =

#TESTS_ENVIRONMENT = @BACKEND_TEST_CMD@
#TESTS_ENVIRONMENT = $(top_srcdir)/test/runtest $(top_builddir)/bin/qhmc

TEST_EXTENSIONS = .lua
LUA_LOG_COMPILER = $(srcdir)/runtest $(top_builddir)/bin/qhmc

TESTS = \
 #fieldIO.lua \
 norm2.lua \
 point.lua \
 subset.lua \
 #wilsonHMC.lua

EXTRA_DIST = \
 fieldIO.lua \
 norm2.lua \
 point.lua \
 subset.lua \
 wilsonHMC.lua \
 fieldIO.ref \
 norm2.ref \
 point.ref \
 subset.ref

distclean-local:
	-rm -f *.out
