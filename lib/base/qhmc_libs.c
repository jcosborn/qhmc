#include "qhmc_internal.h"

void
qhmc_init_libs(int *argc, char ***argv)
{
#ifdef HAVE_QOPQDP
  qhmc_init_qopqdp(argc, argv);
#endif
}

void
qhmc_fini_libs(void)
{
#ifdef HAVE_QOPQDP
  qhmc_fini_qopqdp();
#endif
}

void
qhmc_open_libs(lua_State* L)
{
  luaopen_lfs(L);
  luaopen_bc(L);
  open_qhmc_complex(L);
#ifdef HAVE_QOPQDP
  qhmc_open_qopqdp(L);
#endif
#if 1  // avoid loading 'lfs' since it is statically linked
  int rc = luaL_dostring(L, "table.insert(package.searchers,1, \
    function(m) local function ldr(x) return _ENV[x] end \
    if(m=='lfs') then return ldr else return nil end end)");
  if(rc!=0) {
    printf("error: %s: luaL_dostring failed\n", __func__);
    exit(1);
  }
#endif
}
