#include "qhmc_common.h"

void
get_bool_array(lua_State *L, int idx, int n, int *a)
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qassert(n==lua_objlen(L, idx));
  for(int i=0; i<n; i++) {
    lua_rawgeti(L, idx, i+1);
    a[i] = lua_toboolean(L, -1);
    lua_pop(L, 1);
  }
}

void
get_int_array(lua_State *L, int idx, int n, int *a)
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qassert(n==lua_objlen(L, idx));
  for(int i=0; i<n; i++) {
    lua_rawgeti(L, idx, i+1);
    a[i] = lua_tointeger(L, -1);
    lua_pop(L, 1);
  }
}

void
push_int_array(lua_State *L, int n, int *a)
{
  lua_createtable(L, n, 0);
  for(int i=0; i<n; i++) {
    lua_pushinteger(L, a[i]);
    lua_rawseti(L, -2, i+1);
  }
}

void
get_float_array(lua_State *L, int idx, int n, float *a)
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qassert(n==lua_objlen(L, idx));
  for(int i=0; i<n; i++) {
    lua_rawgeti(L, idx, i+1);
    a[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
}

void
get_double_array(lua_State *L, int idx, int n, double *a)
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qassert(n==lua_objlen(L, idx));
  for(int i=0; i<n; i++) {
    lua_rawgeti(L, idx, i+1);
    a[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
}

void
push_double_array(lua_State *L, int n, double *a)
{
  lua_createtable(L, n, 0);
  for(int i=0; i<n; i++) {
    lua_pushnumber(L, a[i]);
    lua_rawseti(L, -2, i+1);
  }
}

void
push_complex_array(lua_State *L, int n, qhmc_complex_t *a)
{
  lua_createtable(L, n, 0);
  for(int i=0; i<n; i++) {
    qhmc_complex_create(L, a[i].r, a[i].i);
    lua_rawseti(L, -2, i+1);
  }
}

const char *
tableGetString(lua_State *L, int idx, char *key)
{
  lua_getfield(L, idx, key);
  const char *s = luaL_checkstring(L, -1);
  lua_pop(L, 1);
  return s;
}
