#include "qhmc_common.h"

const char *
qhmc_opt_string(lua_State *L, int *idx, int required, char *def)
{
  const char *x = def;
  if(required) {
    lua_pushvalue(L, *idx);
    x = luaL_checkstring(L, *idx);
    lua_pop(L, 1);
    (*idx)++;
  } else {
    if(lua_isstring(L, *idx)) {
      lua_pushvalue(L, *idx);
      x = lua_tostring(L, *idx);
      lua_pop(L, 1);
      (*idx)++;
    }
  }
  return x;
}

int
qhmc_opt_int(lua_State *L, int *idx, int required, int def)
{
  int x = def;
  if(required) {
    x = luaL_checkint(L, *idx);
    (*idx)++;
  } else {
    if(lua_isnumber(L, *idx)) {
      x = lua_tointeger(L, *idx);
      (*idx)++;
    }
  }
  return x;
}

double
qhmc_opt_double(lua_State *L, int *idx, int required, double def)
{
  double x = def;
  if(required) {
    x = luaL_checknumber(L, *idx);
    (*idx)++;
  } else {
    if(lua_isnumber(L, *idx)) {
      x = lua_tonumber(L, *idx);
      (*idx)++;
    }
  }
  return x;
}

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
push_float_array(lua_State *L, int n, float *a)
{
  lua_createtable(L, n, 0);
  for(int i=0; i<n; i++) {
    lua_pushnumber(L, a[i]);
    lua_rawseti(L, -2, i+1);
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
