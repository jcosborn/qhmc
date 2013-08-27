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

// 1: wquark
// 2: coords table
// 3: color
// 4: spin
// 5: re part
// 6: im part
static int
qopqdp_wquark_point(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==6);
  wquark_t *q = qopqdp_wquark_check(L, 1);
  int nd; get_table_len(L, 2, &nd);
  int point[nd]; get_int_array(L, 2, nd, point);
  int color = luaL_checkint(L, 3);
  int spin = luaL_checkint(L, 4);
  double re = luaL_checknumber(L, 5);
  double im = luaL_checknumber(L, 6);
  int node = QDP_node_number(point);
  if(node==QDP_this_node) {
    int index = QDP_index(point);
    QLA_DiracFermion *qdf = QDP_site_ptr_readwrite_D(q->df, index);
    QLA_c_eq_r_plus_ir(QLA_elem_D(*qdf,color,spin), re, im);
  }
  return 0;
}

// 1: wquark
// 2: (optional) subset
static int
qopqdp_wquark_norm2(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  wquark_t *q = qopqdp_wquark_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    const char *s = luaL_checkstring(L, 2);
    if(!strcmp(s,"timeslices")) sub = NULL;
    else sub = qopqdp_check_subset(L, 2);
  }
  if(sub) {
    QLA_Real nrm2;
    QDP_r_eq_norm2_D(&nrm2, q->df, sub);
    lua_pushnumber(L, nrm2);
  } else { // timeslices
    int nt = QDP_coord_size(3);
    QLA_Real nrm2[nt];
    QDP_Subset *ts = qhmcqdp_get_timeslices();
    QDP_r_eq_norm2_D_multi(nrm2, q->df, ts, nt);
    push_double_array(L, nt, nrm2);
  }
  return 1;
}

// 1: wquark
// 2: wquark
// 3: (optional) subset
static int
qopqdp_wquark_redot(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==2 || narg==3);
  wquark_t *q1 = qopqdp_wquark_check(L, 1);
  wquark_t *q2 = qopqdp_wquark_check(L, 2);
  QDP_Subset sub = QDP_all;
  if(narg!=2) {
    const char *s = luaL_checkstring(L, 3);
    if(!strcmp(s,"timeslices")) sub = NULL;
    else sub = qopqdp_check_subset(L, 3);
  }
  if(sub) {
    QLA_Real redot;
    QDP_r_eq_re_D_dot_D(&redot, q1->df, q2->df, sub);
    lua_pushnumber(L, redot);
  } else { // timeslices
    int nt = QDP_coord_size(3);
    QLA_Real redot[nt];
    QDP_Subset *ts = qhmcqdp_get_timeslices();
    QDP_r_eq_re_D_dot_D_multi(redot, q1->df, q2->df, ts, nt);
    push_double_array(L, nt, redot);
  }
  return 1;
}

// 1: wquark
// 2: wquark
// 3: (optional) subset
static int
qopqdp_wquark_imdot(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==2 || narg==3);
  wquark_t *q1 = qopqdp_wquark_check(L, 1);
  wquark_t *q2 = qopqdp_wquark_check(L, 2);
  QDP_Subset sub = QDP_all;
  if(narg!=2) {
    const char *s = luaL_checkstring(L, 3);
    if(!strcmp(s,"timeslices")) sub = NULL;
    else sub = qopqdp_check_subset(L, 3);
  }
  if(sub) {
    QLA_Complex dot;
    QDP_c_eq_D_dot_D(&dot, q1->df, q2->df, sub);
    lua_pushnumber(L, QLA_imag(dot));
  } else { // timeslices
    int nt = QDP_coord_size(3);
    QLA_Complex dot[nt];
    QDP_Subset *ts = qhmcqdp_get_timeslices();
    QDP_c_eq_D_dot_D_multi(dot, q1->df, q2->df, ts, nt);
    double imdot[nt];
    for(int i=0; i<nt; i++) imdot[i] = QLA_imag(dot[i]);
    push_double_array(L, nt, imdot);
  }
  return 1;
}

// 1: wquark
// 2: wquark
// 3: (optional) subset
static int
qopqdp_wquark_dot(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==2 || narg==3);
  wquark_t *q1 = qopqdp_wquark_check(L, 1);
  wquark_t *q2 = qopqdp_wquark_check(L, 2);
  QDP_Subset sub = QDP_all;
  if(narg!=2) {
    const char *s = luaL_checkstring(L, 3);
    if(!strcmp(s,"timeslices")) sub = NULL;
    else sub = qopqdp_check_subset(L, 3);
  }
  if(sub) {
    QLA_Complex dot;
    QDP_c_eq_D_dot_D(&dot, q1->df, q2->df, sub);
    qhmc_complex_create(L, QLA_real(dot), QLA_imag(dot));
  } else { // timeslices
    int nt = QDP_coord_size(3);
    QLA_Complex dot[nt];
    QDP_Subset *ts = qhmcqdp_get_timeslices();
    QDP_c_eq_D_dot_D_multi(dot, q1->df, q2->df, ts, nt);
    qhmc_complex_t qdot[nt];
    for(int i=0; i<nt; i++) {
      qdot[i].r = QLA_real(dot[i]);
      qdot[i].i = QLA_imag(dot[i]);
    }
    push_complex_array(L, nt, qdot);
  }
  return 1;
}

static int
qopqdp_wquark_combine(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==3 || narg==4);
  wquark_t *qd = qopqdp_wquark_check(L, 1); //output quark field
  int nqs; get_table_len(L, 2, &nqs); // input quark fields
  wquark_t *qs[nqs]; qopqdp_wquark_array_check(L, 2, nqs, qs);
  int nc; get_table_len(L, 3, &nc); // input coefficients
  qassert(nqs==nc);
  double c[nc]; get_double_array(L, 3, nc, c);
  QLA_Real qc[nc]; for(int i=0; i<nc; i++) qc[i] = c[i];
  QDP_Subset sub = QDP_all;
  if(narg>3) {
    sub = qopqdp_check_subset(L, 4);
  }
  QDP_D_eq_r_times_D(qd->df, &qc[0], qs[0]->df, sub);
  for(int i=1; i<nqs; i++) {
    QDP_D_peq_r_times_D(qd->df, &qc[i], qs[i]->df, sub);
  }
  return 0;
}

// 1: wquark
// 2: gamma
// 3: wquark
// TODO 4: (optional) options ("peq")
static int
qopqdp_wquark_gamma(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==3);
  wquark_t *qd = qopqdp_wquark_check(L, 1); //output quark field
  int g = luaL_checkint(L, 2);
  wquark_t *qs = qopqdp_wquark_check(L, 3); //output quark field
  QDP_Subset sub = QDP_all;
  QDP_D_eq_gamma_times_D(qd->df, qs->df, g, sub);
  return 0;
}

// 1: wquark
// 2: gauge
// 3: r
// 4: (optional) n
static int
qopqdp_wquark_smearGauss(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==3 || narg==4);
  wquark_t *q = qopqdp_wquark_check(L, 1);
  gauge_t *g = qopqdp_gauge_check(L, 2);
  double r = luaL_checknumber(L, 3);
  int n = luaL_optint(L, 4, 1);

  QDP_DiracFermion *t, *tt[3], *ts[3], *tf[3], *tb1[3], *tb2[3];
  QDP_ShiftDir fwd[3], bck[3];
  t = QDP_create_D();
  for(int i=0; i<3; i++) {
    tt[i] = t;
    ts[i] = q->df;
    tf[i] = QDP_create_D();
    tb1[i] = QDP_create_D();
    tb2[i] = QDP_create_D();
    fwd[i] = QDP_forward;
    bck[i] = QDP_backward;
  }
  QLA_Real s0 = 1.0/(1.0+6*r);
  QLA_Real s1 = s0*r;
  for(int i=0; i<n; i++) {
    QDP_D_eq_r_times_D(t, &s1, q->df, QDP_all);
    QDP_D_veq_sD(tf, tt, QDP_neighbor, fwd, QDP_all, 3);
    QDP_D_veq_Ma_times_D(tb1, g->links, tt, QDP_all, 3);
    QDP_D_veq_sD(tb2, tb1, QDP_neighbor, bck, QDP_all, 3);
    QDP_D_eq_r_times_D(q->df, &s0, q->df, QDP_all);
    QDP_D_vpeq_M_times_D(ts, g->links, tf, QDP_all, 3);
    QDP_D_vpeq_D(ts, tb2, QDP_all, 3);
    for(int i=0; i<3; i++) {
      QDP_discard_D(tf[i]);
      QDP_discard_D(tb2[i]);
    }
  }
  QDP_destroy_D(t);
  for(int i=0; i<3; i++) {
    QDP_destroy_D(tf[i]);
    QDP_destroy_D(tb1[i]);
    QDP_destroy_D(tb2[i]);
  }
  return 0;
}

static struct luaL_Reg wquark_reg[] = {
  { "__gc",       qopqdp_wquark_gc },
  { "zero",       qopqdp_wquark_zero },
  { "random",     qopqdp_wquark_random },
  { "randomU1",   qopqdp_wquark_randomU1 },
  { "set",        qopqdp_wquark_set },
  { "point",      qopqdp_wquark_point },
  { "norm2",      qopqdp_wquark_norm2 },
  { "Re_dot",     qopqdp_wquark_redot },
  { "Im_dot",     qopqdp_wquark_imdot },
  { "dot",        qopqdp_wquark_dot },
  { "combine",    qopqdp_wquark_combine },
  { "gamma",      qopqdp_wquark_gamma },
  { "smearGauss", qopqdp_wquark_smearGauss },
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
