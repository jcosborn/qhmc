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

#define STR(x) STR2(x)
#define STR2(x) # x

#define if0 if(1)
#define printf0(...) do { if0 printf(__VA_ARGS__); } while(0)
#define printerr(...) fprintf(stderr, __VA_ARGS__)
#define printerr0(...) do { if0 printerr(__VA_ARGS__); } while(0)
#define TRACE do{if0{printf("%s %s %i\n", __FILE__, __func__, __LINE__);fflush(stdout);}} while(0)
#define ABORT(c) exit(c)
#define LABORT(L,c) luaL_error(L, "aborted with code %d\n", c)
#define qerror(c,...)    do {TRACE;printf(__VA_ARGS__); ABORT(c);  } while(0)
#define qlerror(L,c,...) do {TRACE;printf(__VA_ARGS__);LABORT(L,c);} while(0)
#define qassert(x)    do{if(!(x)) qerror(1,"assert failed (%s)\n",#x);}while(0)
#define qlassert(L,x) do{if(!(x))qlerror(L,1,"assert failed (%s)\n",#x);}while(0)
#define qerror0(c,...)    do { if0  qerror(c,__VA_ARGS__); } while(0)
#define qlerror0(L,c,...) do { if0 qlerror(L,c,__VA_ARGS__); } while(0)
#define qassert0(x) do { if0 qassert(x); } while(0)
#define qlassert0(L,x) do { if0 qlassert(L,x); } while(0)

#define get_table_len(L, idx, n) luaL_checktype(L, idx, LUA_TTABLE); *(n) = lua_objlen(L, idx)
#define tableLength(L, idx) lua_rawlen(L, idx)
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

#define BEGIN_ARGS int nextarg=1, nargs=lua_gettop(L)
#define GET_STRING(v) const char *v=qhmc_opt_string(L,&nextarg,1,NULL)
#define OPT_STRING(v,d) const char *v=qhmc_opt_string(L,&nextarg,0,d)
#define GET_INT(v) int v=qhmc_opt_int(L,&nextarg,1,0)
#define OPT_INT(v,d) int v=qhmc_opt_int(L,&nextarg,0,d)
#define GET_INT_ARRAY(n,v) int n; get_table_len(L,nextarg,&n); int v[n]; qhmc_get_int_array(L,nextarg,n,v); nextarg++
#define GET_DOUBLE(v) double v=qhmc_opt_double(L,&nextarg,1,0)
#define OPT_DOUBLE(v,d) double v=qhmc_opt_double(L,&nextarg,0,d)
#define GET_DOUBLE_ARRAY(n,v) int n; get_table_len(L,nextarg,&n); double v[n]; qhmc_get_double_array(L,nextarg,n,v); nextarg++
#define OPT_FUNCTION_INDEX(f,d) int f=(d);if(lua_isfunction(L,nextarg))f=nextarg++
#define GET_TABLE_LEN_INDEX(l,i) int l,i=nextarg; get_table_len(L,i,&l); nextarg++
#define OPT_TABLE_LEN_INDEX(l,i) int l=-1,i=0; { if(lua_istable(L,nextarg)) { i=nextarg; get_table_len(L,i,&l); nextarg++; } }
#define END_ARGS qlassert(L,nextarg==nargs+1)
#define GET_AS_DOUBLE_ARRAY(n,t) int n=qhmc_opt_as_double_array_len(L,nextarg,1,0); double t[n==0?1:abs(n)]; qhmc_opt_as_double_array(L,&nextarg,1,n,t,0,NULL)
#define OPT_AS_DOUBLE_ARRAY(n,t,dn,dt) int n=qhmc_opt_as_double_array_len(L,nextarg,0,dn); double t[n==0?1:abs(n)]; qhmc_opt_as_double_array(L,&nextarg,0,n,t,dn,dt)

#define QHMC_NEWUSERDATA(t, v) \
  t *v = lua_newuserdata(L, sizeof(t)); \
  QHMC_USERTABLE_CREATE(L, -1); \
  QHMC_PTRTABLE_SET(L, -1, v);

#define QHMC_USERTABLE_CREATE(L,i) { int _i=lua_absindex(L,i); lua_newtable(L); lua_setuservalue(L,_i); }
#define QHMC_USERTABLE_GETFIELD(L,i,k) lua_getuservalue(L,i); lua_getfield(L,-1,k); lua_remove(L,-2)
#define QHMC_USERTABLE_SETFIELD(L,i,k) lua_getuservalue(L,i); lua_insert(L,-2); lua_setfield(L,-2,k); lua_pop(L,1)

#define QHMC_PTRTABLE_SET(L, i, p) {			\
	int _i = lua_absindex(L, i);			\
	lua_getfield(L, LUA_REGISTRYINDEX, "ptrtable"); \
	lua_pushvalue(L, _i);				\
	lua_rawsetp(L, -2, p);				\
	lua_pop(L, 1);					\
      }

#define QHMC_PTRTABLE_GET(L, p) {			\
    lua_getfield(L, LUA_REGISTRYINDEX, "ptrtable");	\
    lua_rawgetp(L, -1, p);				\
    lua_remove(L, -2);					\
  }


const char *qhmc_opt_string(lua_State *L, int *idx, int required, char *def);
int qhmc_opt_int(lua_State *L, int *idx, int required, int def);
double qhmc_opt_double(lua_State *L, int *idx, int required, double def);
int qhmc_opt_as_double_array_len(lua_State *L, int idx, int required, int def);
void qhmc_opt_as_double_array(lua_State *L, int *idx, int required, int n, double *t, int dn, double *def);

void qhmc_get_bool_array(lua_State *L, int idx, int n, int *a);
void qhmc_get_int_array(lua_State *L, int idx, int n, int *a);
void qhmc_push_int_array(lua_State *L, int n, int *a);
void qhmc_get_float_array(lua_State *L, int idx, int n, float *a);
void qhmc_push_float_array(lua_State *L, int n, float *a);
void qhmc_get_double_array(lua_State *L, int idx, int n, double *a);
void qhmc_push_double_array(lua_State *L, int n, double *a);
void qhmc_push_complex_array(lua_State *L, int n, qhmc_complex_t *a);
const char *qhmc_tableGetString(lua_State *L, int idx, char *key);
#define tableGetString qhmc_tableGetString

void qhmc_open_qhmc(lua_State *L);
