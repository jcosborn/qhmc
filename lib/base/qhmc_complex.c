/// Complex number type.
// The constructor (complex) is in the global namespace (@{global.complex}).
//
// Complex numbers support basic operations +,-,*,/ among themselves or with
// regular numbers.
//
// The real or imaginary part can be accessed with .r or .i (i.e. z.r or z.i).
// This works for setting and getting the value.
// @classmod complex
#include "qhmc_internal.h"
#include <math.h>

static char *mtname = "complex";

qhmc_complex_t *
qhmc_complex_check(lua_State *L, int idx)
{
  qhmc_complex_t *c = luaL_checkudata(L, idx, mtname);
  return c;
}

qhmc_complex_t *
qhmc_complex_test(lua_State *L, int idx)
{
  return luaL_testudata(L, idx, mtname);
}

qhmc_complex_t *
qhmc_opt_complex(lua_State *L, int *idx, int required, qhmc_complex_t *def)
{
  qhmc_complex_t *t;
  if(required) {
    t = luaL_checkudata(L, *idx, mtname);
    (*idx)++;
  } else {
    t = luaL_testudata(L, *idx, mtname);
    if(t==NULL) t = def;
    else (*idx)++;
  }
  return t;
}

qhmc_complex_t *
qhmc_opt_as_complex_ptr(lua_State *L, int *idx, int required,
			qhmc_complex_t *buf, qhmc_complex_t *def)
{
  qhmc_complex_t *t;
  if(lua_isnumber(L, *idx)) {
    t = buf;
    t->r = lua_tonumber(L, *idx);
    t->i = 0;
    (*idx)++;
  } else {
    if(required) {
      t = luaL_checkudata(L, *idx, mtname);
      (*idx)++;
    } else {
      t = luaL_testudata(L, *idx, mtname);
      if(t==NULL) t = def;
      else (*idx)++;
    }
  }
  return t;
}

void
qhmc_opt_as_complex_def_real(lua_State *L, int *idx, int required,
			     qhmc_complex_t *buf, double def)
{
  if(lua_isnumber(L, *idx)) {
    buf->r = lua_tonumber(L, *idx);
    buf->i = 0;
    (*idx)++;
  } else {
    if(required) {
      qhmc_complex_t *t = luaL_checkudata(L, *idx, mtname);
      buf->r = t->r;
      buf->i = t->i;
      (*idx)++;
    } else {
      qhmc_complex_t *t = luaL_testudata(L, *idx, mtname);
      if(t==NULL) {
	buf->r = def;
	buf->i = 0;
      } else {
	buf->r = t->r;
	buf->i = t->i;
	(*idx)++;
      }
    }
  }
}

void
qhmc_complex_get_as(lua_State *L, int idx, qhmc_complex_t *c)
{
  if(lua_isnumber(L,idx)) {
    c->r = lua_tonumber(L, idx);
    c->i = 0;
  } else {
    qhmc_complex_t *z = qhmc_complex_check(L, idx);
    c->r = z->r;
    c->i = z->i;
  }
}

// returns len:
//  -1 if single value found
//   n if table of 'n' values found
// def otherwise
int
qhmc_opt_as_complex_array_len(lua_State *L, int idx, int required, int def)
{
  int len = def, tlen = -1, i = idx;
  int type = lua_type(L, idx);
  if(type==LUA_TTABLE) {
    tlen = tableLength(L, idx);
    tableGetIndex(L, idx, 1);
    i = -1;
  }
  qhmc_complex_t *t, b;
  int ii = i;
  t = qhmc_opt_as_complex_ptr(L, &ii, required, &b, NULL);
  if(t!=NULL) len = tlen;
  if(type==LUA_TTABLE) {
    lua_pop(L, 1);
  }
  return len;
}

void
qhmc_opt_as_complex_array(lua_State *L, int *idx, int required, int n,
			  qhmc_complex_t *t, int dn, qhmc_complex_t *def)
{
  int type = lua_type(L, *idx);
  if(type==LUA_TTABLE) {
    tableGetIndex(L, *idx, 1);
    qhmc_complex_t *tt, b;
    int ii = -1;
    tt = qhmc_opt_as_complex_ptr(L, &ii, required, &b, NULL);
    lua_pop(L, 1);
    if(tt==NULL) {
      for(int i=0; i<abs(n); i++) t[i] = def[i];
    } else {
      for(int i=0; i<n; i++) {
        tableGetIndex(L, *idx, i+1);
	qhmc_complex_get_as(L, -1, &t[i]);
        lua_pop(L, 1);
      }
      (*idx)++;
    }
  } else {
    qhmc_opt_as_complex_ptr(L, idx, required, t, def);
  }
}

static int
qhmc_complex_add(lua_State *L)
{
  qhmc_complex_t *c = qhmc_complex_create(L, 0, 0);
  if(lua_isnumber(L,1)) { // second must be complex
    qhmc_complex_t *c2 = qhmc_complex_check(L,2);
    c->r = c2->r + lua_tonumber(L,1);
    c->i = c2->i;
  } else {
    qhmc_complex_t *c1 = qhmc_complex_check(L,1);
    if(lua_isnumber(L,2)) {
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
  if(lua_isnumber(L,1)) { // second must be complex
    qhmc_complex_t *c2 = qhmc_complex_check(L,2);
    c->r = lua_tonumber(L,1) - c2->r;
    c->i = - c2->i;
  } else {
    qhmc_complex_t *c1 = qhmc_complex_check(L,1);
    if(lua_isnumber(L,2)) {
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
  if(lua_isnumber(L,1)) { // second must be complex
    qhmc_complex_t *c = qhmc_complex_create(L, 0, 0);
    qhmc_complex_t *c2 = qhmc_complex_check(L,2);
    double r = lua_tonumber(L,1);
    c->r = r*c2->r;
    c->i = r*c2->i;
  } else {
    qhmc_complex_t *c1 = qhmc_complex_check(L,1);
    if(lua_isnumber(L,2)) {
      qhmc_complex_t *c = qhmc_complex_create(L, 0, 0);
      double r = lua_tonumber(L,2);
      c->r = r*c1->r;
      c->i = r*c1->i;
    } else {
      qhmc_complex_t *c2 = qhmc_complex_test(L,2);
      if(c2==NULL) { // unknown object, use it's metamethod
	luaL_getmetafield(L, 2, "__mul");
	lua_pushvalue(L, 2);
	lua_pushvalue(L, 1);
	lua_call(L, 2, 1);
      } else {
	qhmc_complex_t *c = qhmc_complex_create(L, 0, 0);
	c->r = c1->r*c2->r - c1->i*c2->i;
	c->i = c1->r*c2->i + c1->i*c2->r;
      }
    }
  }
  return 1;
}

static int
qhmc_complex_div(lua_State *L)
{
  qhmc_complex_t *c = qhmc_complex_create(L, 0, 0);
  if(lua_isnumber(L,1)) { // second must be complex
    qhmc_complex_t *c2 = qhmc_complex_check(L,2);
    double r = lua_tonumber(L,1);
    double s = r/(c2->r*c2->r + c2->i*c2->i);
    c->r = s*c2->r;
    c->i = -s*c2->i;
  } else {
    qhmc_complex_t *c1 = qhmc_complex_check(L,1);
    if(lua_isnumber(L,2)) {
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
qhmc_complex_eq(lua_State *L)
{
  qhmc_complex_t t1, t2;
  qhmc_complex_t *z1 = qhmc_opt_as_complex_ptr(L, (int[]){1}, 0, &t1, NULL);
  qhmc_complex_t *z2 = qhmc_opt_as_complex_ptr(L, (int[]){2}, 0, &t2, NULL);
  printf("eq: %p  %p\n", z1, z2);
  if(z1==NULL || z2==NULL) {
    lua_pushboolean(L, 0);
  } else {
    lua_pushboolean(L, (z1->r==z2->r)&&(z1->i==z2->i));
  }
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
  qlerror(L, 1, "invalid complex element (%s)\n", key);
  return 0;
}

// need to add format options
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

/// Set complex number from another, or 1 or 2 real numbers.
//  @function set
//  @tparam ?number|complex x real part or complex value
//  @tparam[opt=0] number y imaginary part (only if x is not complex)
static int
qhmc_complex_set(lua_State *L)
{
  qhmc_complex_t *c = qhmc_complex_check(L, 1);
  qhmc_complex_get_as(L, 2, c);
  return 0;
}

/// Return complex conjugate.
//  @function conj
//  @treturn complex conjugate of self
static int
qhmc_complex_conj(lua_State *L)
{
  qhmc_complex_t *c = qhmc_complex_check(L, 1);
  qhmc_complex_create(L, c->r, -c->i);
  return 1;
}

/// Add another number to this one (+=).
//  @function peq
//  @tparam ?real|complex z number to add to self
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

/// Subtract another number from this one (-=).
//  @function meq
//  @tparam ?real|complex z number to subtract from self
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

/// Return square root.
//  @function sqrt
//  @treturn complex square root of self
static int
qhmc_complex_sqrt(lua_State *L)
{
  qhmc_complex_t *c = qhmc_complex_check(L, 1);
  double r = c->r;
  double i = c->i;
  double n = r*r + i*i;
  double sr = sqrt(fabs(0.5*(n+r)));
  double si = copysign(sqrt(fabs(0.5*(n-r))), i);
  qhmc_complex_create(L, sr, si);
  return 1;
}

/// Return absolute value squared.
//  @function norm2
//  @treturn complex absolute value squared of self
static int
qhmc_complex_norm2(lua_State *L)
{
  qhmc_complex_t *c = qhmc_complex_check(L, 1);
  double r = c->r;
  double i = c->i;
  double n = r*r + i*i;
  lua_pushnumber(L, n);
  return 1;
}

static struct luaL_Reg qhmc_complex_mt_reg[] = {
  { "__add",      qhmc_complex_add },
  { "__sub",      qhmc_complex_sub },
  { "__mul",      qhmc_complex_mul },
  { "__div",      qhmc_complex_div },
  //{ "__pow",      qhmc_complex_pow },
  { "__unm",      qhmc_complex_unm },
  { "__eq",       qhmc_complex_eq },
  { "__index",    qhmc_complex_index },
  { "__newindex", qhmc_complex_newindex },
  { "__tostring", qhmc_complex_tostring },
  { "set",        qhmc_complex_set },
  { "conj",       qhmc_complex_conj },
  { "peq",        qhmc_complex_peq },
  { "meq",        qhmc_complex_meq },
  { "sqrt",       qhmc_complex_sqrt },
  { "norm2",      qhmc_complex_norm2 },
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

// doc in global.cdoc
/// Create complex number.
//  @function complex
//  @tparam ?number|complex x real part or complex value
//  @tparam[opt=0] number y imaginary part (only if x is not complex)
//  @treturn complex
static int
qhmc_complex(lua_State* L)
{
  int nargs = lua_gettop(L);
  qlassert(L, nargs<=2);
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
