#ifdef lua_c
#include "qhmc.h"
#define main lmain
#define luaL_openlibs(L) luaL_openlibs(L); qhmc_open_libs(L)
#endif

int lmain(int argc, char **argv);
