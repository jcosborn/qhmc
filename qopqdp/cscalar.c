#include <string.h>
#include <math.h>
#include "qhmc_qopqdp_common.h"

static char *mtname = "qopqdp.cscalar";

cscalar_t *
qopqdp_cscalar_check(lua_State *L, int idx)
{
  luaL_checkudata(L, idx, mtname);
  cscalar_t *cs = lua_touserdata(L, idx);
  return cs;
}

void
qopqdp_cscalar_array_check(lua_State *L, int idx, int n, cscalar_t *s[n])
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qassert(lua_objlen(L, idx)==n);
  lua_pushvalue(L, idx); // make copy for indexing convenience
  for(int i=0; i<n; i++) {
    lua_pushnumber(L, i+1);
    lua_gettable(L, -2);
    s[i] = qopqdp_cscalar_check(L, -1);
    lua_pop(L, 1);
  }
  lua_pop(L, 1);
}

static void
qopqdp_cscalar_free(lua_State *L, int idx)
{
  cscalar_t *s = qopqdp_cscalar_check(L, idx);
  QDP_destroy_C(s->c);
}

static int
qopqdp_cscalar_gc(lua_State *L)
{
  qopqdp_cscalar_free(L, -1);
  return 0;
}

static int
qopqdp_cscalar_zero(lua_State *L)
{
  BEGIN_ARGS;
  GET_CSCALAR(s);
  OPT_SUBSET(sub, s->lat, QDP_all_L(s->qlat));
  END_ARGS;
  QDP_C_eq_zero(s->c, sub);
  return 0;
}

static int
qopqdp_cscalar_point(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==5);
  cscalar_t *s = qopqdp_cscalar_check(L, 1);
  int nd; get_table_len(L, 2, &nd);
  int point[nd]; get_int_array(L, 2, nd, point);
  double re = luaL_checknumber(L, 3);
  double im = luaL_checknumber(L, 4);
  int node = QDP_node_number(point);
  if(node==QDP_this_node) {
    int index = QDP_index(point);
    QLA_Complex *qc = QDP_site_ptr_readwrite_C(s->c, index);
    QLA_c_eq_r_plus_ir(*qc, re, im);
  }
  return 0;
}

static int
qopqdp_cscalar_random(lua_State *L)
{
  BEGIN_ARGS;
  GET_CSCALAR(s);
  OPT_SUBSET(sub, s->lat, QDP_all_L(s->qlat));
  END_ARGS;
  QDP_C_eq_gaussian_S(s->c, qopqdp_srs, sub);
  QLA_Real r = sqrt(0.5); // normalize to sigma^2 = 1/2
  QDP_C_eq_r_times_C(s->c, &r, s->c, sub);
  return 0;
}

static void
lnormalize_C(QLA_Complex *x, int i)
{
  QLA_Real n = QLA_norm2_c(*x);
  if(n!=0) n = 1/sqrt(n);
  QLA_c_eq_r_times_c(*x, n, *x);
}

static int
qopqdp_cscalar_randomU1(lua_State *L)
{
  BEGIN_ARGS;
  GET_CSCALAR(s);
  OPT_SUBSET(sub, s->lat, QDP_all_L(s->qlat));
  END_ARGS;
  QDP_C_eq_gaussian_S(s->c, qopqdp_srs, sub);
  QDP_C_eq_funci(s->c, lnormalize_C, sub);
  return 0;
}

static int
qopqdp_cscalar_set(lua_State *L)
{
  BEGIN_ARGS;
  GET_CSCALAR(s1);
  GET_CSCALAR(s2);
  OPT_SUBSET(sub, s1->lat, QDP_all_L(s1->qlat));
  END_ARGS;
  QDP_C_eq_C(s1->c, s2->c, sub);
  return 0;
}

static int
qopqdp_cscalar_combine(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==3 || narg==4);
  cscalar_t *sd = qopqdp_cscalar_check(L, 1);
  int nss; get_table_len(L, 2, &nss);
  cscalar_t *ss[nss]; qopqdp_cscalar_array_check(L, 2, nss, ss);
  int nc; get_table_len(L, 3, &nc);
  qassert(nss==nc);
  double c[nc]; get_double_array(L, 3, nc, c);
  QLA_Real qc[nc]; for(int i=0; i<nc; i++) qc[i] = c[i];
  lattice_t *lat = sd->lat;
  QDP_Subset sub = QDP_all_L(lat->qlat);
  if(narg>3) {
    sub = qopqdp_check_subset(L, 4, lat);
  }
  QDP_C_eq_r_times_C(sd->c, &qc[0], ss[0]->c, sub);
  for(int i=1; i<nss; i++) {
    QDP_C_peq_r_times_C(sd->c, &qc[i], ss[i]->c, sub);
  }
  return 0;
}

static int
qopqdp_cscalar_norm2(lua_State *L)
{
  BEGIN_ARGS;
  GET_CSCALAR(s);
  OPT_SUBSETS(subs, ns, s->lat, QDP_all_and_empty_L(s->qlat), 1);
  END_ARGS;
  if(ns==1) {
    QLA_Real nrm2;
    QDP_r_eq_norm2_C(&nrm2, s->c, subs[0]);
    lua_pushnumber(L, nrm2);
  } else {
    QLA_Real nrm2[ns];
    QDP_r_eq_norm2_C_multi(nrm2, s->c, subs, ns);
    push_real_array(L, ns, nrm2);
  }
  return 1;
}

static int
qopqdp_cscalar_redot(lua_State *L)
{
  BEGIN_ARGS;
  GET_CSCALAR(s1);
  GET_CSCALAR(s2);
  OPT_SUBSETS(subs, ns, s1->lat, QDP_all_and_empty_L(s1->qlat), 1);
  END_ARGS;
  if(ns==1) {
    QLA_Real redot;
    QDP_r_eq_re_C_dot_C(&redot, s1->c, s2->c, subs[0]);
    lua_pushnumber(L, redot);
  } else {
    QLA_Real redot[ns];
    QDP_r_eq_re_C_dot_C_multi(redot, s1->c, s2->c, subs, ns);
    push_real_array(L, ns, redot);
  }
  return 1;
}

static struct luaL_Reg cscalar_reg[] = {
  { "__gc",    qopqdp_cscalar_gc },
  { "zero",    qopqdp_cscalar_zero },
  { "point",   qopqdp_cscalar_point },
  { "random",  qopqdp_cscalar_random },
  { "randomU1",qopqdp_cscalar_randomU1 },
  { "set",     qopqdp_cscalar_set },
  { "combine", qopqdp_cscalar_combine },
  { "norm2",   qopqdp_cscalar_norm2 },
  { "Re_dot",  qopqdp_cscalar_redot },
  { NULL, NULL}
};

cscalar_t *
qopqdp_cscalar_create(lua_State* L, lattice_t *lat)
{
  if(lat==NULL) lat = qopqdp_get_default_lattice(L);
  cscalar_t *s = lua_newuserdata(L, sizeof(cscalar_t));
  s->c = QDP_create_C_L(lat->qlat);
  s->lat = lat;
  s->qlat = lat->qlat;
  QDP_C_eq_zero(s->c, QDP_all_L(lat->qlat));
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, cscalar_reg);
  }
  lua_setmetatable(L, -2);
  return s;
}
