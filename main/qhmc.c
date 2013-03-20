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
open_qhmc(lua_State* L) {
  luaopen_lfs(L);
  open_qopqdp(L);
  int rc = luaL_dostring(L, "table.insert(package.searchers,1, \
    function(m) return function(m) return _ENV[m] end end)");
  if(rc!=0) {
    printf("error: %s: luaL_dostring failed\n", __func__);
    exit(1);
  }
}

int
main(int argc, char *argv[])
{
  init_libs(&argc, &argv);

  lmain(argc, argv);

  fini_libs();

  return 0;
}
