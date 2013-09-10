#include <stdlib.h>
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
#include "qhmc_complex.h"

#if 0
#define lua_objlen lua_rawlen
#define luaL_register(L,n,l) do { \
    if(n) { luaL_newlibtable(L,l); lua_pushvalue(L,-1); lua_setglobal(L,n); } \
    luaL_setfuncs(L,l,0); \
  } while(0)
#endif

#define if0 if(1)
#define printf0(...) do { if0 printf(__VA_ARGS__); } while(0)
#define printerr(...) fprintf(stderr, __VA_ARGS__)
#define printerr0(...) do { if0 printerr(__VA_ARGS__); } while(0)
#define TRACE printf0("%s %s %i\n", __FILE__, __func__, __LINE__)
#define ABORT(code) exit(code)
#define LABORT(L, code) luaL_error(L, "aborted with code %i\n", code)
#define qerror(...) do { TRACE;	printf(__VA_ARGS__); ABORT(1); } while(0)
#define qlerror(L,...) do{ TRACE; printf(__VA_ARGS__); LABORT(L, 1); }while(0)
#define qassert(x) do{ if(!(x)) qerror("assert failed (%s)\n", #x); }while(0)
#define qlassert(L,x) do{if(!(x))qlerror(L,"assert failed (%s)\n",#x);}while(0)
#define qerror0(...) do { if0 qerror(__VA_ARGS__); } while(0)
#define qlerror0(L, ...) do { if0 qlerror(L, __VA_ARGS__); } while(0)
#define qassert0(x) do { if0 qassert(x); } while(0)
#define qlassert0(L, x) do { if0 qlassert(L, x); } while(0)

#define get_table_len(L, idx, n) luaL_checktype(L, idx, LUA_TTABLE); *(n) = lua_objlen(L, idx)
#define tableLoopKeys(L, idx)		\
  qassert(lua_type(L,idx)==LUA_TTABLE); \
  lua_pushnil(L);			\
  while(lua_next(L,idx)) {		\
    lua_pushvalue(L, -2);		\
    const char *tableLoopKey = lua_tostring(L, -1);
#define tableLoopKeysIfKeySetInt(L, key, var)			\
    if(strcmp(tableLoopKey,key)==0) (var) = lua_tointeger(L,-2)
#define tableLoopKeysIfKeySetDouble(L, key, var)			\
    if(strcmp(tableLoopKey,key)==0) (var) = lua_tonumber(L,-2)
#define tableLoopKeysEnd(L) \
  /*else {					\
      /-* error *-/				\
    }*/ \
    lua_pop(L, 2); \
  }
#define tableGetField(L, idx, key) lua_getfield(L, idx, key)
#define tableGetIndex(L, idx, key) lua_rawgeti(L, idx, key)

#define START_ARGS int nextarg=1, nargs=lua_gettop(L)
#define GET_DOUBLE(v) double v=lua_tonumber(L,nextarg); nextarg++
#define GET_TABLE_LEN_INDEX(l,i) int l,i=nextarg; get_table_len(L,i,&l); nextarg++
#define END_ARGS qassert(nextarg==nargs+1)

void get_bool_array(lua_State *L, int idx, int n, int *a);
void get_int_array(lua_State *L, int idx, int n, int *a);
void push_int_array(lua_State *L, int n, int *a);
void get_float_array(lua_State *L, int idx, int n, float *a);
void get_double_array(lua_State *L, int idx, int n, double *a);
void push_double_array(lua_State *L, int n, double *a);
void push_complex_array(lua_State *L, int n, qhmc_complex_t *a);
const char *tableGetString(lua_State *L, int idx, char *key);
