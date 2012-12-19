#include <string.h>
#include <math.h>
#include "qhmc_qopqdp_common.h"

static char *mtname = "qopqdp.dwquark";

dwquark_t *
qopqdp_dwquark_check(lua_State *L, int idx)
{
  luaL_checkudata(L, idx, mtname);
  dwquark_t *q = lua_touserdata(L, idx);
  return q;
}

void
qopqdp_dwquark_array_check(lua_State *L, int idx, int n, dwquark_t *q[n])
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qassert(lua_objlen(L, idx)==n);
  lua_pushvalue(L, idx); // make copy for indexing convenience
  for(int i=0; i<n; i++) {
    lua_pushnumber(L, i+1);
    lua_gettable(L, -2);
    q[i] = qopqdp_dwquark_check(L, -1);
    lua_pop(L, 1);
  }
  lua_pop(L, 1);
}

static void
qopqdp_dwquark_free(lua_State *L, int idx)
{
  dwquark_t *q = qopqdp_dwquark_check(L, idx);
  for(int s=0; s<q->ls; s++) {
    QDP_destroy_D(q->df[s]);
  }
  free(q->df);
}

static int
qopqdp_dwquark_gc(lua_State *L)
{
  qopqdp_dwquark_free(L, -1);
  return 0;
}

static int
qopqdp_dwquark_zero(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  dwquark_t *q = qopqdp_dwquark_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    sub = qopqdp_check_subset(L, 2);
  }
  for(int s=0; s<q->ls; s++) {
    QDP_D_eq_zero(q->df[s], sub);
  }
  return 0;
}

static int
qopqdp_dwquark_random(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  dwquark_t *q = qopqdp_dwquark_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    sub = qopqdp_check_subset(L, 2);
  }
  QLA_Real r = sqrt(0.5); // normalize to sigma^2 = 1/2
  for(int s=0; s<q->ls; s++) {
    QDP_D_eq_gaussian_S(q->df[s], qopqdp_srs, sub);
    QDP_D_eq_r_times_D(q->df[s], &r, q->df[s], sub);
  }
  return 0;
}

static void
lnormalize_D(NCPROT QLA_DiracFermion(*x), int i)
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
qopqdp_dwquark_randomU1(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  dwquark_t *q = qopqdp_dwquark_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    sub = qopqdp_check_subset(L, 2);
  }
  for(int s=0; s<q->ls; s++) {
    QDP_D_eq_gaussian_S(q->df[s], qopqdp_srs, sub);
    QDP_D_eq_funci(q->df[s], lnormalize_D, sub);
  }
  return 0;
}

static int
qopqdp_dwquark_set(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==2 || narg==3);
  dwquark_t *q1 = qopqdp_dwquark_check(L, 1);
  dwquark_t *q2 = qopqdp_dwquark_check(L, 2);
  qassert(q1->ls==q2->ls);
  QDP_Subset sub = QDP_all;
  if(narg!=2) {
    sub = qopqdp_check_subset(L, 3);
  }
  for(int s=0; s<q1->ls; s++) {
    QDP_D_eq_D(q1->df[s], q2->df[s], sub);
  }
  return 0;
}

static int
qopqdp_dwquark_norm2(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  dwquark_t *q = qopqdp_dwquark_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    sub = qopqdp_check_subset(L, 2);
  }
  QLA_Real nrm2=0;
  for(int s=0; s<q->ls; s++) {
    QLA_Real tnrm2;
    QDP_r_eq_norm2_D(&tnrm2, q->df[s], sub);
    nrm2 += tnrm2;
  }
  lua_pushnumber(L, nrm2);
  return 1;
}

static int
qopqdp_dwquark_redot(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==2 || narg==3);
  dwquark_t *q1 = qopqdp_dwquark_check(L, 1);
  dwquark_t *q2 = qopqdp_dwquark_check(L, 2);
  qassert(q1->ls==q2->ls);
  QDP_Subset sub = QDP_all;
  if(narg!=2) {
    sub = qopqdp_check_subset(L, 3);
  }
  QLA_Real redot=0;
  for(int s=0; s<q1->ls; s++) {
    QLA_Real tredot;
    QDP_r_eq_re_D_dot_D(&tredot, q1->df[s], q2->df[s], sub);
    redot += tredot;
  }
  lua_pushnumber(L, redot);
  return 1;
}

static int
qopqdp_dwquark_combine(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==3 || narg==4);
  dwquark_t *qd = qopqdp_dwquark_check(L, 1);
  int nqs; get_table_len(L, 2, &nqs);
  dwquark_t *qs[nqs]; qopqdp_dwquark_array_check(L, 2, nqs, qs);
  for(int i=0; i<nqs; i++) qassert(qd->ls==qs[i]->ls);
  int nc; get_table_len(L, 3, &nc);
  qassert(nqs==nc);
  double c[nc]; get_double_array(L, 3, nc, c);
  QLA_Real qc[nc]; for(int i=0; i<nc; i++) qc[i] = c[i];
  QDP_Subset sub = QDP_all;
  if(narg>3) {
    sub = qopqdp_check_subset(L, 4);
  }
  for(int s=0; s<qd->ls; s++) {
    QDP_D_eq_r_times_D(qd->df[s], &qc[0], qs[0]->df[s], sub);
    for(int i=1; i<nqs; i++) {
      QDP_D_peq_r_times_D(qd->df[s], &qc[i], qs[i]->df[s], sub);
    }
  }
  return 0;
}

static struct luaL_Reg dwquark_reg[] = {
  { "__gc",     qopqdp_dwquark_gc },
  { "zero",     qopqdp_dwquark_zero },
  { "random",   qopqdp_dwquark_random },
  { "randomU1", qopqdp_dwquark_randomU1 },
  { "set",      qopqdp_dwquark_set },
  { "norm2",    qopqdp_dwquark_norm2 },
  { "Re_dot",   qopqdp_dwquark_redot },
  { "combine",  qopqdp_dwquark_combine },
  { NULL, NULL}
};

dwquark_t *
qopqdp_dwquark_create(lua_State* L, int ls)
{
  dwquark_t *q = lua_newuserdata(L, sizeof(dwquark_t));
  q->ls = ls;
  q->df = malloc(ls*sizeof(*q->df));
  for(int s=0; s<ls; s++) {
    q->df[s] = QDP_create_D();
    QDP_D_eq_zero(q->df[s], QDP_all);
  }
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, dwquark_reg);
  }
  lua_setmetatable(L, -2);
  return q;
}
