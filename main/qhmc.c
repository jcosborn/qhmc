#include <stdlib.h>
#include <stdio.h>
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
#include "qhmc.h"
#include "qhmclibs.h"

void
init_libs(int *argc, char ***argv)
{
  init_qopqdp(argc, argv);
}

void
fini_libs(void)
{
  fini_qopqdp();
}

void
open_qhmc(lua_State* L)
{
  luaopen_lfs(L);
  open_qhmc_complex(L);
  open_qopqdp(L);
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

int
main(int argc, char *argv[])
{
  init_libs(&argc, &argv);

  lmain(argc, argv);

  fini_libs();

  return 0;
}
