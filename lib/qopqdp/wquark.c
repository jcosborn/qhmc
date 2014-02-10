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
qopqdp_wquark_clone(lua_State *L)
{
  BEGIN_ARGS;
  GET_WQUARK(q2);
  OPT_SUBSET(sub, q2->lat, QDP_all_L(q2->qlat));
  END_ARGS;
  wquark_t *q1 = qopqdp_wquark_create_unset(L, q2->nc, q2->lat);
  QDP_D_eq_D(q1->df, q2->df, sub);
  return 1;
}

static int
qopqdp_wquark_zero(lua_State *L)
{
  BEGIN_ARGS;
  GET_WQUARK(q);
  OPT_SUBSET(sub, q->lat, QDP_all_L(q->qlat));
  END_ARGS;
  QDP_D_eq_zero(q->df, sub);
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

static int
qopqdp_wquark_random(lua_State *L)
{
  BEGIN_ARGS;
  GET_WQUARK(q);
  OPT_SUBSET(sub, q->lat, QDP_all_L(q->qlat));
  END_ARGS;
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
  BEGIN_ARGS;
  GET_WQUARK(q);
  OPT_SUBSET(sub, q->lat, QDP_all_L(q->qlat));
  END_ARGS;
  QDP_D_eq_gaussian_S(q->df, qopqdp_srs, sub);
  QDP_D_eq_funci(q->df, lnormalize_D, sub);
  return 0;
}

static int
qopqdp_wquark_set(lua_State *L)
{
  BEGIN_ARGS;
  GET_WQUARK(q1);
  GET_WQUARK(q2);
  OPT_SUBSET(sub, q1->lat, QDP_all_L(q1->qlat));
  END_ARGS;
  QDP_D_eq_D(q1->df, q2->df, sub);
  return 0;
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
  //double c[nc]; get_double_array(L, 3, nc, c);
  //QLA_Real qc[nc]; for(int i=0; i<nc; i++) qc[i] = c[i];
  QDP_Subset sub = QDP_all;
  if(narg>3) {
    sub = qopqdp_check_subset(L, 4, qd->lat);
  }
  for(int i=0; i<nqs; i++) {
    lua_rawgeti(L, 3, i+1);
    if(lua_type(L,-1)==LUA_TNUMBER) {
      QLA_Real c = lua_tonumber(L, -1);
      if(i==0) {
	QDP_D_eq_r_times_D(qd->df, &c, qs[0]->df, sub);
      } else {
	QDP_D_peq_r_times_D(qd->df, &c, qs[i]->df, sub);
      }
    } else { // assume complex
      qhmc_complex_t *qc = qhmc_complex_check(L, -1);
      QLA_Complex c;
      QLA_c_eq_r_plus_ir(c, qc->r, qc->i);
      if(i==0) {
	QDP_D_eq_c_times_D(qd->df, &c, qs[0]->df, sub);
      } else {
	QDP_D_peq_c_times_D(qd->df, &c, qs[i]->df, sub);
      }
    }
  }
  return 0;
}

static int
qopqdp_wquark_getSite(lua_State *L)
{
#define NC QDP_get_nc(q->df)
  BEGIN_ARGS;
  GET_WQUARK(q);
  GET_TABLE_LEN_INDEX(nd,ip);
  END_ARGS;
  int site[nd]; get_int_array(L, ip, nd, site);
  QLA_DiracFermion qdf;
  int node = QDP_node_number_L(q->qlat, site);
  if(node==QDP_this_node) {
    int index = QDP_index_L(q->qlat, site);
    QLA_DiracFermion *qdfp = QDP_site_ptr_readonly_D(q->df, index);
    QLA_D_eq_D(&qdf, qdfp);
  } else {
    QLA_D_eq_zero(&qdf);
  }
  QMP_sum_double_array((double *)&qdf, 4*2*QLA_Nc);
  push_complex_array(L, 4*QLA_Nc, (qhmc_complex_t *)&qdf);
  return 1;
#undef NC
}

// 1: wquark
// 2: (optional) subset
static int
qopqdp_wquark_norm2(lua_State *L)
{
  BEGIN_ARGS;
  GET_WQUARK(q);
  OPT_SUBSETS(subs, ns, q->lat, QDP_all_and_empty_L(q->qlat), 1);
  END_ARGS;
  if(ns==1) {
    QLA_Real nrm2;
    QDP_r_eq_norm2_D(&nrm2, q->df, subs[0]);
    lua_pushnumber(L, nrm2);
  } else {
    QLA_Real nrm2[ns];
    QDP_r_eq_norm2_D_multi(nrm2, q->df, subs, ns);
    push_double_array(L, ns, nrm2);
  }
  return 1;
}

// 1: wquark
// 2: wquark
// 3: (optional) subset
static int
qopqdp_wquark_redot(lua_State *L)
{
  BEGIN_ARGS;
  GET_WQUARK(q1);
  GET_WQUARK(q2);
  OPT_SUBSETS(subs, ns, q1->lat, QDP_all_and_empty_L(q1->qlat), 1);
  END_ARGS;
  if(ns==1) {
    QLA_Real redot;
    QDP_r_eq_re_D_dot_D(&redot, q1->df, q2->df, subs[0]);
    lua_pushnumber(L, redot);
  } else {
    QLA_Real redot[ns];
    QDP_r_eq_re_D_dot_D_multi(redot, q1->df, q2->df, subs, ns);
    push_double_array(L, ns, redot);
  }
  return 1;
}

// 1: wquark
// 2: wquark
// 3: (optional) subset
static int
qopqdp_wquark_imdot(lua_State *L)
{
  BEGIN_ARGS;
  GET_WQUARK(q1);
  GET_WQUARK(q2);
  OPT_SUBSETS(subs, ns, q1->lat, QDP_all_and_empty_L(q1->qlat), 1);
  END_ARGS;
  if(ns==1) {
    QLA_Complex dot;
    QDP_c_eq_D_dot_D(&dot, q1->df, q2->df, subs[0]);
    lua_pushnumber(L, QLA_imag(dot));
  } else {
    QLA_Complex dot[ns];
    QDP_c_eq_D_dot_D_multi(dot, q1->df, q2->df, subs, ns);
    double imdot[ns];
    for(int i=0; i<ns; i++) imdot[i] = QLA_imag(dot[i]);
    push_double_array(L, ns, imdot);
  }
  return 1;
}

// 1: wquark
// 2: wquark
// 3: (optional) subset
static int
qopqdp_wquark_dot(lua_State *L)
{
  BEGIN_ARGS;
  GET_WQUARK(q1);
  GET_WQUARK(q2);
  OPT_SUBSETS(subs, ns, q1->lat, QDP_all_and_empty_L(q1->qlat), 1);
  END_ARGS;
  if(ns==1) {
    QLA_Complex dot;
    QDP_c_eq_D_dot_D(&dot, q1->df, q2->df, subs[0]);
    qhmc_complex_create(L, QLA_real(dot), QLA_imag(dot));
  } else {
    QLA_Complex dot[ns];
    QDP_c_eq_D_dot_D_multi(dot, q1->df, q2->df, subs, ns);
    qhmc_complex_t qdot[ns];
    for(int i=0; i<ns; i++) {
      qdot[i].r = QLA_real(dot[i]);
      qdot[i].i = QLA_imag(dot[i]);
    }
    push_complex_array(L, ns, qdot);
  }
  return 1;
}

// 1: wquark
// 2: gamma
// 3: wquark
// 4: (optional) subset
// TODO 5: (optional) options ("peq")
static int
qopqdp_wquark_gamma(lua_State *L)
{
  BEGIN_ARGS;
  GET_WQUARK(qd);
  GET_INT(g);
  GET_WQUARK(qs);
  OPT_SUBSET(sub, qd->lat, QDP_all_L(qd->qlat));
  END_ARGS;
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
#define NC QDP_get_nc(q->df)
  BEGIN_ARGS;
  GET_WQUARK(q);
  GET_GAUGE(g);
  GET_DOUBLE(r);
  OPT_INT(n,1);
  END_ARGS;
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
#undef NC
}

static struct luaL_Reg wquark_reg[] = {
  { "__gc",       qopqdp_wquark_gc },
  { "clone",      qopqdp_wquark_clone },
  { "zero",       qopqdp_wquark_zero },
  { "point",      qopqdp_wquark_point },
  { "random",     qopqdp_wquark_random },
  { "randomU1",   qopqdp_wquark_randomU1 },
  { "set",        qopqdp_wquark_set },
  { "combine",    qopqdp_wquark_combine },
  { "getSite",    qopqdp_wquark_getSite },
  { "norm2",      qopqdp_wquark_norm2 },
  { "Re_dot",     qopqdp_wquark_redot },
  { "Im_dot",     qopqdp_wquark_imdot },
  { "dot",        qopqdp_wquark_dot },
  { "gamma",      qopqdp_wquark_gamma },
  { "smearGauss", qopqdp_wquark_smearGauss },
  { NULL, NULL}
};

wquark_t *
qopqdp_wquark_create_unset(lua_State *L, int nc, lattice_t *lat)
{
#define NC nc
  if(nc==0) nc = QOPQDP_DEFAULTNC;
  if(lat==NULL) lat = qopqdp_get_default_lattice(L);
  wquark_t *q = lua_newuserdata(L, sizeof(wquark_t));
  q->df = QDP_create_D_L(lat->qlat);
  q->lat = lat;
  q->qlat = lat->qlat;
  q->nc = nc;
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, wquark_reg);
  }
  lua_setmetatable(L, -2);
  return q;
#undef NC
}

wquark_t *
qopqdp_wquark_create(lua_State *L, int nc, lattice_t *lat)
{
  wquark_t *q = qopqdp_wquark_create_unset(L, nc, lat);
  QDP_D_eq_zero(q->df, QDP_all);
  return q;
}
