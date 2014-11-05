#include "qhmc_internal.h"
#include <string.h>
#include <libgen.h>

static char **g_argv = NULL;

void
qhmc_init_libs(int *argc, char ***argv)
{
#ifdef HAVE_QOPQDP
  qhmc_init_qopqdp(argc, argv);
#endif
#ifdef HAVE_QUDA
  qhmc_init_quda(argc, argv);
#endif
  g_argv = *argv;
}

void
qhmc_fini_libs(void)
{
#ifdef HAVE_QOPQDP
  qhmc_fini_qopqdp();
#endif
#ifdef HAVE_QUDA
  qhmc_fini_quda();
#endif
}

static void
addPath(lua_State* L, char *s)
{
  const int size = 1024;
  char buf[size];
  int rc = snprintf(buf,size, "package.path = '%s/?.lua;' .. package.path", s);
  if(rc<0||rc>=size) {
    printf("error %s:%s:%i: snprintf failed\n", __FILE__, __func__, __LINE__);
    exit(1);
  }
  //printf("%s\n", buf);
  rc = luaL_dostring(L, buf);
  if(rc!=0) {
    printf("error %s:%s:%i: luaL_dostring failed\n", __FILE__, __func__, __LINE__);
    exit(1);
  }
}

void
qhmc_open_libs(lua_State* L)
{
  qhmc_open_qhmc(L);
  luaopen_lfs(L);
  luaopen_bc(L);
  open_qhmc_complex(L);
#ifdef HAVE_QOPQDP
  qhmc_open_qopqdp(L);
#endif
#ifdef HAVE_QUDA
  qhmc_open_quda(L);
#endif

  { // add search paths for libraries
    // srcdir
    addPath(L, SRCDIRLUA);
    // installdir
    addPath(L, PREFIXLUA);
    // exedir
    char *dirc;
    dirc = strdup(g_argv[0]);
    addPath(L, dirname(dirc));
    free(dirc);
    // scriptdir?
    // other
    addPath(L, "./lua");
    addPath(L, ".");
  }

  { // avoid loading 'lfs' since it is statically linked
    int rc = luaL_dostring(L, "table.insert(package.searchers,1, \
      function(m) local function ldr(x) return _ENV[x] end	 \
      if(m=='lfs') then return ldr else return nil end end)");
    if(rc!=0) {
      printf("error %s:%s:%i: luaL_dostring failed\n", __FILE__, __func__, __LINE__);
      exit(1);
    }
  }

  if(qhmc_master()) {
    if(0) { // print package search path
      int rc = luaL_dostring(L, "print(\"package.path = \"..package.path)");
      if(rc!=0) {
	printf("error %s:%s:%i: luaL_dostring failed\n", __FILE__, __func__, __LINE__);
	exit(1);
      }
    }
    { // print path to loaded packages
      int rc = luaL_dostring(L, "local ps; ps=function(m)\
 for k,f in ipairs(package.searchers) do if f~=ps then r,e=f(m);	\
 if type(r)==\"function\" then print(\"#loading: \"..e);return r,e end	\
 end end return nil end; table.insert(package.searchers,1,ps)");
      if(rc!=0) {
	printf("error %s:%s:%i: luaL_dostring failed\n", __FILE__, __func__, __LINE__);
	exit(1);
      }
    }
  }

  { // load Init.lua
    int rc = luaL_dostring(L, "require 'Init'");
    if(rc!=0) {
      printf("error %s:%s:%i: luaL_dostring failed\n", __FILE__, __func__, __LINE__);
      exit(1);
    }
  }

}
