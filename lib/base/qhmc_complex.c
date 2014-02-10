#include "qhmc_internal.h"
#include <math.h>

static char *mtname = "complex";

qhmc_complex_t *
qhmc_complex_check(lua_State *L, int idx)
{
  qhmc_complex_t *c = luaL_checkudata(L, idx, mtname);
  return c;
}

void
qhmc_complex_get_as(lua_State *L, int idx, qhmc_complex_t *c)
{
  int type = lua_type(L, idx);
  if(type==LUA_TNUMBER) {
    c->r = lua_tonumber(L, idx);
    c->i = 0;
  } else {
    qhmc_complex_t *z = qhmc_complex_check(L, idx);
    c->r = z->r;
    c->i = z->i;
  }
}

static int
qhmc_complex_add(lua_State *L)
{
  qhmc_complex_t *c = qhmc_complex_create(L, 0, 0);
  if(lua_type(L,1)==LUA_TNUMBER) { // second must be complex
    qhmc_complex_t *c2 = qhmc_complex_check(L,2);
    c->r = c2->r + lua_tonumber(L,1);
    c->i = c2->i;
  } else {
    qhmc_complex_t *c1 = qhmc_complex_check(L,1);
    if(lua_type(L,2)==LUA_TNUMBER) {
      c->r = c1->r + lua_tonumber(L,2);
      c->i = c1->i;
    } else {
      qhmc_complex_t *c2 = qhmc_complex_check(L,2);
      c->r = c1->r + c2->r;
      c->i = c1->i + c2->i;
    }
  }
  return 1;
}

static int
qhmc_complex_sub(lua_State *L)
{
  qhmc_complex_t *c = qhmc_complex_create(L, 0, 0);
  if(lua_type(L,1)==LUA_TNUMBER) { // second must be complex
    qhmc_complex_t *c2 = qhmc_complex_check(L,2);
    c->r = lua_tonumber(L,1) - c2->r;
    c->i = - c2->i;
  } else {
    qhmc_complex_t *c1 = qhmc_complex_check(L,1);
    if(lua_type(L,2)==LUA_TNUMBER) {
      c->r = c1->r - lua_tonumber(L,2);
      c->i = c1->i;
    } else {
      qhmc_complex_t *c2 = qhmc_complex_check(L,2);
      c->r = c1->r - c2->r;
      c->i = c1->i - c2->i;
    }
  }
  return 1;
}

static int
qhmc_complex_mul(lua_State *L)
{
  qhmc_complex_t *c = qhmc_complex_create(L, 0, 0);
  if(lua_type(L,1)==LUA_TNUMBER) { // second must be complex
    qhmc_complex_t *c2 = qhmc_complex_check(L,2);
    double r = lua_tonumber(L,1);
    c->r = r*c2->r;
    c->i = r*c2->i;
  } else {
    qhmc_complex_t *c1 = qhmc_complex_check(L,1);
    if(lua_type(L,2)==LUA_TNUMBER) {
      double r = lua_tonumber(L,2);
      c->r = r*c1->r;
      c->i = r*c1->i;
    } else {
      qhmc_complex_t *c2 = qhmc_complex_check(L,2);
      c->r = c1->r*c2->r - c1->i*c2->i;
      c->i = c1->r*c2->i + c1->i*c2->r;
    }
  }
  return 1;
}

static int
qhmc_complex_div(lua_State *L)
{
  qhmc_complex_t *c = qhmc_complex_create(L, 0, 0);
  if(lua_type(L,1)==LUA_TNUMBER) { // second must be complex
    qhmc_complex_t *c2 = qhmc_complex_check(L,2);
    double r = lua_tonumber(L,1);
    double s = r/(c2->r*c2->r + c2->i*c2->i);
    c->r = s*c2->r;
    c->i = -s*c2->i;
  } else {
    qhmc_complex_t *c1 = qhmc_complex_check(L,1);
    if(lua_type(L,2)==LUA_TNUMBER) {
      double r = lua_tonumber(L,2);
      double s = 1/r;
      c->r = s*c1->r;
      c->i = s*c1->i;
    } else {
      qhmc_complex_t *c2 = qhmc_complex_check(L,2);
      double s = 1/(c2->r*c2->r + c2->i*c2->i);
      c->r = s*(c1->r*c2->r + c1->i*c2->i);
      c->i = s*(c1->i*c2->r - c1->r*c2->i);
    }
  }
  return 1;
}

static int
qhmc_complex_unm(lua_State *L)
{
  qhmc_complex_t *c1 = qhmc_complex_check(L,1);
  qhmc_complex_create(L, -c1->r, -c1->i);
  return 1;
}

static int
qhmc_complex_index(lua_State *L)
{
  size_t len;
  const char *key = lua_tolstring(L, 2, &len);
  if(len==1) {
    if(key[0]=='r') {
      qhmc_complex_t *c = qhmc_complex_check(L, 1);
      lua_pushnumber(L, c->r);
      return 1;
    } else if(key[0]=='i') {
      qhmc_complex_t *c = qhmc_complex_check(L, 1);
      lua_pushnumber(L, c->i);
      return 1;
    }
  }
  lua_getmetatable(L, 1);
  lua_getfield(L, -1, key);
  return 1;
}

static int
qhmc_complex_newindex(lua_State *L)
{
  size_t len;
  const char *key = lua_tolstring(L, 2, &len);
  if(len==1) {
    if(key[0]=='r') {
      qhmc_complex_t *c = qhmc_complex_check(L, 1);
      c->r = lua_tonumber(L, 3);
      return 0;
    } else if(key[0]=='i') {
      qhmc_complex_t *c = qhmc_complex_check(L, 1);
      c->i = lua_tonumber(L, 3);
      return 0;
    }
  }
  qlerror0(L,1,"invalid complex element (%s)\n", key);
  return 0;
}

static int
qhmc_complex_tostring(lua_State *L)
{
  qhmc_complex_t *c = qhmc_complex_check(L, 1);
  lua_pushfstring(L, "(%f,%f)", (lua_Number)c->r, (lua_Number)c->i);
  //if(c->i<0) {
  //  lua_pushfstring(L, "%f%fi", (lua_Number)c->r, (lua_Number)c->i);
  //} else {
  //  lua_pushfstring(L, "%f+%fi", (lua_Number)c->r, (lua_Number)c->i);
  //}
  return 1;
}

static int
qhmc_complex_set(lua_State *L)
{
  qhmc_complex_t *c = qhmc_complex_check(L, 1);
  qhmc_complex_get_as(L, 2, c);
  return 0;
}

static int
qhmc_complex_conj(lua_State *L)
{
  qhmc_complex_t *c = qhmc_complex_check(L, 1);
  qhmc_complex_create(L, c->r, -c->i);
  return 1;
}

static int
qhmc_complex_peq(lua_State *L)
{
  qhmc_complex_t c2;
  qhmc_complex_get_as(L, 2, &c2);
  qhmc_complex_t *c = qhmc_complex_check(L, 1);
  c->r += c2.r;
  c->i += c2.i;
  return 0;
}

static int
qhmc_complex_meq(lua_State *L)
{
  qhmc_complex_t c2;
  qhmc_complex_get_as(L, 2, &c2);
  qhmc_complex_t *c = qhmc_complex_check(L, 1);
  c->r -= c2.r;
  c->i -= c2.i;
  return 0;
}

static struct luaL_Reg qhmc_complex_mt_reg[] = {
  { "__add",      qhmc_complex_add },
  { "__sub",      qhmc_complex_sub },
  { "__mul",      qhmc_complex_mul },
  { "__div",      qhmc_complex_div },
  //{ "__pow",      qhmc_complex_pow },
  { "__unm",      qhmc_complex_unm },
  //{ "__eq",       qhmc_complex_eq },
  { "__index",    qhmc_complex_index },
  { "__newindex", qhmc_complex_newindex },
  { "__tostring", qhmc_complex_tostring },
  { "set",        qhmc_complex_set },
  { "conj",       qhmc_complex_conj },
  { "peq",        qhmc_complex_peq },
  { "meq",        qhmc_complex_meq },
  { NULL, NULL}
};

qhmc_complex_t *
qhmc_complex_create(lua_State* L, double re, double im)
{
  qhmc_complex_t *c = lua_newuserdata(L, sizeof(qhmc_complex_t));
  c->r = re;
  c->i = im;
  if(luaL_newmetatable(L, mtname)) {
    //lua_pushvalue(L, -1);
    //lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, qhmc_complex_mt_reg);
  }
  lua_setmetatable(L, -2);
  return c;
}

static int
qhmc_complex(lua_State* L)
{
  int nargs = lua_gettop(L);
  qassert(nargs<=2);
  qhmc_complex_t *c = qhmc_complex_create(L, 0, 0);
  if(nargs==1) {
    qhmc_complex_get_as(L, 1, c);
  } else { // nargs==2
    c->r = lua_tonumber(L, 1);
    c->i = lua_tonumber(L, 2);
  }
  return 1;
}

void
open_qhmc_complex(lua_State* L)
{
  lua_pushcfunction(L, qhmc_complex);
  lua_setglobal(L, "complex");
}
