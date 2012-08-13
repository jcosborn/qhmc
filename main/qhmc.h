#ifdef lua_c
#define main lmain
#define luaL_openlibs(L) luaL_openlibs(L); open_qhmc(L)
#endif

extern void open_qhmc(lua_State* L);
extern int lmain(int argc, char **argv);
