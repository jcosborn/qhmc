#include <string.h>
#include <math.h>
#include "qhmc_qopqdp_common.h"

static char *mtname = "qopqdp.wquark";

wquark_t *
qopqdp_wquark_check(lua_State *L, int idx)
{
  luaL_checkudata(L, idx, mtname);
  wquark_t *q = lua_touserdata(L, idx);
  return q;
}

void
qopqdp_wquark_array_check(lua_State *L, int idx, int n, wquark_t *q[n])
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qassert(lua_objlen(L, idx)==n);
  lua_pushvalue(L, idx); // make copy for indexing convenience
  for(int i=0; i<n; i++) {
    lua_pushnumber(L, i+1);
    lua_gettable(L, -2);
    q[i] = qopqdp_wquark_check(L, -1);
    lua_pop(L, 1);
  }
  lua_pop(L, 1);
}

static void
qopqdp_wquark_free(lua_State *L, int idx)
{
  wquark_t *q = qopqdp_wquark_check(L, idx);
  QDP_destroy_D(q->df);
}

static int
qopqdp_wquark_gc(lua_State *L)
{
  qopqdp_wquark_free(L, -1);
  return 0;
}

static int
qopqdp_wquark_zero(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  wquark_t *q = qopqdp_wquark_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    sub = qopqdp_check_subset(L, 2);
  }
  QDP_D_eq_zero(q->df, sub);
  return 0;
}

static int
qopqdp_wquark_random(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  wquark_t *q = qopqdp_wquark_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    sub = qopqdp_check_subset(L, 2);
  }
  QDP_D_eq_gaussian_S(q->df, qopqdp_srs, sub);
  QLA_Real r = sqrt(0.5); // normalize to sigma^2 = 1/2
  QDP_D_eq_r_times_D(q->df, &r, q->df, sub);
  return 0;
}

static void
lnormalize_D(QLA_DiracFermion *x, int i)
{
  for(int ic=0; ic<QLA_Nc; ic++) {
    for(int is=0; is<QLA_Ns; is++) {
      QLA_Complex z = QLA_elem_D(*x,ic,is);
      QLA_Real n = QLA_norm2_c(z);
      if(n!=0) n = 1/sqrt(n);
      QLA_c_eq_r_times_c(QLA_elem_D(*x,ic,is), n, z);
    }
  }
}

static int
qopqdp_wquark_randomU1(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  wquark_t *q = qopqdp_wquark_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    sub = qopqdp_check_subset(L, 2);
  }
  QDP_D_eq_gaussian_S(q->df, qopqdp_srs, sub);
  QDP_D_eq_funci(q->df, lnormalize_D, sub);
  return 0;
}

static int
qopqdp_wquark_set(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==2 || narg==3);
  wquark_t *q1 = qopqdp_wquark_check(L, 1);
  wquark_t *q2 = qopqdp_wquark_check(L, 2);
  QDP_Subset sub = QDP_all;
  if(narg!=2) {
    sub = qopqdp_check_subset(L, 3);
  }
  QDP_D_eq_D(q1->df, q2->df, sub);
  return 0;
}

static int
qopqdp_wquark_norm2(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  wquark_t *q = qopqdp_wquark_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    sub = qopqdp_check_subset(L, 2);
  }
  QLA_Real nrm2;
  QDP_r_eq_norm2_D(&nrm2, q->df, sub);
  lua_pushnumber(L, nrm2);
  return 1;
}

static int
qopqdp_wquark_redot(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==2 || narg==3);
  wquark_t *q1 = qopqdp_wquark_check(L, 1);
  wquark_t *q2 = qopqdp_wquark_check(L, 2);
  QDP_Subset sub = QDP_all;
  if(narg!=2) {
    sub = qopqdp_check_subset(L, 3);
  }
  QLA_Real redot;
  QDP_r_eq_re_D_dot_D(&redot, q1->df, q2->df, sub);
  lua_pushnumber(L, redot);
  return 1;
}

static int
qopqdp_wquark_combine(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==3 || narg==4);
  wquark_t *qd = qopqdp_wquark_check(L, 1);
  int nqs; get_table_len(L, 2, &nqs);
  wquark_t *qs[nqs]; qopqdp_wquark_array_check(L, 2, nqs, qs);
  int nc; get_table_len(L, 3, &nc);
  qassert(nqs==nc);
  double c[nc]; get_double_array(L, 3, nc, c);
  QLA_Real qc[nc]; for(int i=0; i<nc; i++) qc[i] = c[i];
  QDP_Subset sub = QDP_all;
  if(narg>32) {
    sub = qopqdp_check_subset(L, 4);
  }
  QDP_D_eq_r_times_D(qd->df, &qc[0], qs[0]->df, sub);
  for(int i=1; i<nqs; i++) {
    QDP_D_peq_r_times_D(qd->df, &qc[i], qs[i]->df, sub);
  }
  return 0;
}

static struct luaL_Reg wquark_reg[] = {
  { "__gc",     qopqdp_wquark_gc },
  { "zero",     qopqdp_wquark_zero },
  { "random",   qopqdp_wquark_random },
  { "randomU1", qopqdp_wquark_randomU1 },
  { "set",      qopqdp_wquark_set },
  { "norm2",    qopqdp_wquark_norm2 },
  { "Re_dot",   qopqdp_wquark_redot },
  { "combine",  qopqdp_wquark_combine },
  { NULL, NULL}
};

wquark_t *
qopqdp_wquark_create(lua_State* L)
{
  wquark_t *q = lua_newuserdata(L, sizeof(wquark_t));
  q->df = QDP_create_D();
  QDP_D_eq_zero(q->df, QDP_all);
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, wquark_reg);
  }
  lua_setmetatable(L, -2);
  return q;
}
