#include "qhmc_internal.h"
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199309L
#endif
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>

const char *
qhmc_opt_string(lua_State *L, int *idx, int required, char *def)
{
  const char *x = def;
  if(required) {
    lua_pushvalue(L, *idx);
    x = luaL_checkstring(L, -1);
    lua_pop(L, 1);
    (*idx)++;
  } else {
    if(lua_isstring(L, *idx)) {
      lua_pushvalue(L, *idx);
      x = lua_tostring(L, -1);
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

// returns len:
//  -1 if single value found
//   n if table of 'n' values found
// def otherwise
int
qhmc_opt_as_double_array_len(lua_State *L, int idx, int required, int def)
{
  int len = def, tlen = -1, i = idx;
  int type = lua_type(L, idx);
  if(type==LUA_TTABLE) {
    tlen = tableLength(L, idx);
    tableGetIndex(L, idx, 1);
    i = -1;
  }
  int ii = i;
  qhmc_opt_double(L, &ii, required, 0);
  if(ii>i) len = tlen;
  if(type==LUA_TTABLE) {
    lua_pop(L, 1);
  }
  return len;
}

void
qhmc_opt_as_double_array(lua_State *L, int *idx, int required, int n,
			 double *t, int dn, double *def)
{
  int type = lua_type(L, *idx);
  if(type==LUA_TTABLE) {
    tableGetIndex(L, *idx, 1);
    int ii = -1;
    qhmc_opt_double(L, &ii, required, 0);
    lua_pop(L, 1);
    if(ii==-1) {
      for(int i=0; i<abs(n); i++) t[i] = def[i];
    } else {
      for(int i=0; i<n; i++) {
        tableGetIndex(L, *idx, i+1);
	t[i] = lua_tonumber(L, -1);
        lua_pop(L, 1);
      }
      (*idx)++;
    }
  } else {
    double df = 0;
    if(def!=NULL) df = *def;
    *t = qhmc_opt_double(L, idx, required, df);
  }
}

void
qhmc_get_bool_array(lua_State *L, int idx, int n, int *a)
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qlassert(L, n==lua_objlen(L, idx));
  for(int i=0; i<n; i++) {
    lua_rawgeti(L, idx, i+1);
    a[i] = lua_toboolean(L, -1);
    lua_pop(L, 1);
  }
}

void
qhmc_get_int_array(lua_State *L, int idx, int n, int *a)
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qlassert(L, n==lua_objlen(L, idx));
  for(int i=0; i<n; i++) {
    lua_rawgeti(L, idx, i+1);
    a[i] = lua_tointeger(L, -1);
    lua_pop(L, 1);
  }
}

void
qhmc_push_int_array(lua_State *L, int n, int *a)
{
  lua_createtable(L, n, 0);
  for(int i=0; i<n; i++) {
    lua_pushinteger(L, a[i]);
    lua_rawseti(L, -2, i+1);
  }
}

void
qhmc_get_float_array(lua_State *L, int idx, int n, float *a)
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qlassert(L, n==lua_objlen(L, idx));
  for(int i=0; i<n; i++) {
    lua_rawgeti(L, idx, i+1);
    a[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
}

void
qhmc_push_float_array(lua_State *L, int n, float *a)
{
  lua_createtable(L, n, 0);
  for(int i=0; i<n; i++) {
    lua_pushnumber(L, a[i]);
    lua_rawseti(L, -2, i+1);
  }
}

void
qhmc_get_double_array(lua_State *L, int idx, int n, double *a)
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qlassert(L, n==lua_objlen(L, idx));
  for(int i=0; i<n; i++) {
    lua_rawgeti(L, idx, i+1);
    a[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
}

void
qhmc_push_double_array(lua_State *L, int n, double *a)
{
  lua_createtable(L, n, 0);
  for(int i=0; i<n; i++) {
    lua_pushnumber(L, a[i]);
    lua_rawseti(L, -2, i+1);
  }
}

void
qhmc_push_complex_array(lua_State *L, int n, qhmc_complex_t *a)
{
  lua_createtable(L, n, 0);
  for(int i=0; i<n; i++) {
    qhmc_complex_create(L, a[i].r, a[i].i);
    lua_rawseti(L, -2, i+1);
  }
}

const char *
qhmc_tableGetString(lua_State *L, int idx, char *key)
{
  lua_getfield(L, idx, key);
  const char *s = luaL_checkstring(L, -1);
  lua_pop(L, 1);
  return s;
}

static int
qhmc_isMaster(lua_State *L)
{
  qlassert(L, lua_gettop(L)==0);
  int ismaster = qhmc_master();
  lua_pushboolean(L, ismaster);
  return 1;
}

static int
qhmc_dtime(lua_State *L)
{
  qlassert(L, lua_gettop(L)==0);
  double t = 0;
#if _POSIX_TIMERS > 0
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  t = (double) ts.tv_sec + 1e-9 * (double) ts.tv_nsec;
#else
  struct timeval tv;
  gettimeofday(&tv, NULL);
  t = (double) tv.tv_sec + 1e-6 * (double) tv.tv_usec;
#endif
  lua_pushnumber(L, t);
  return 1;
}

static int fdout, fderr;

static int
qhmc_remapout(lua_State *L)
{
  qlassert(L, lua_gettop(L)==1);
  lua_pushvalue(L, -1);
  const char *s = luaL_checkstring(L, -1);
  lua_pop(L, 1);
  int fd = creat(s, 0666);
  fflush(stdout);
  fdout = dup(1);
  dup2(fd, 1);
  fflush(stderr);
  fderr = dup(2);
  dup2(fd, 2);
  close(fd);
  return 0;
}

static int
qhmc_restoreout(lua_State *L)
{
  qlassert(L, lua_gettop(L)==0);
  fflush(stdout);
  dup2(fdout, 1);
  close(fdout);
  fflush(stderr);
  dup2(fderr, 2);
  close(fderr);
  return 0;
}

static struct luaL_Reg qhmc_reg[] = {
  { "isMaster",       qhmc_isMaster },
  { "dtime",          qhmc_dtime },
  { "remapout",       qhmc_remapout },
  { "restoreout",     qhmc_restoreout },
  { NULL, NULL}
};

void
qhmc_open_qhmc(lua_State *L)
{
  luaL_register(L, "qhmc", qhmc_reg);
}
