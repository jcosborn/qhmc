#include <string.h>
#include <math.h>
#include "qhmc_qopqdp_common.h"

static char *mtname = "qopqdp.squark";

squark_t *
qopqdp_squark_check(lua_State *L, int idx)
{
  luaL_checkudata(L, idx, mtname);
  squark_t *q = lua_touserdata(L, idx);
  return q;
}

void
qopqdp_squark_array_check(lua_State *L, int idx, int n, squark_t *q[n])
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qassert(lua_objlen(L, idx)==n);
  lua_pushvalue(L, idx); // make copy for indexing convenience
  for(int i=0; i<n; i++) {
    lua_pushnumber(L, i+1);
    lua_gettable(L, -2);
    q[i] = qopqdp_squark_check(L, -1);
    lua_pop(L, 1);
  }
  lua_pop(L, 1);
}

static void
qopqdp_squark_free(lua_State *L, int idx)
{
  squark_t *q = qopqdp_squark_check(L, idx);
  QDP_destroy_V(q->cv);
}

static int
qopqdp_squark_gc(lua_State *L)
{
  qopqdp_squark_free(L, -1);
  return 0;
}

static int
qopqdp_squark_clone(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  squark_t *q1 = qopqdp_squark_create_unset(L);
  squark_t *q2 = qopqdp_squark_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    sub = qopqdp_check_subset(L, 2);
  }
  QDP_V_eq_V(q1->cv, q2->cv, sub);
  return 1;
}

static int
qopqdp_squark_zero(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  squark_t *q = qopqdp_squark_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    sub = qopqdp_check_subset(L, 2);
  }
  QDP_V_eq_zero(q->cv, sub);
  return 0;
}

static int
qopqdp_squark_getSite(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==2);
  squark_t *q = qopqdp_squark_check(L, 1);
  int nd; get_table_len(L, 2, &nd);
  int site[nd]; get_int_array(L, 2, nd, site);
  QLA_ColorVector qcv;
  int node = QDP_node_number(site);
  if(node==QDP_this_node) {
    int index = QDP_index(site);
    QLA_ColorVector *qcvp = QDP_site_ptr_readonly_V(q->cv, index);
    QLA_V_eq_V(&qcv, qcvp);
  } else {
    QLA_V_eq_zero(&qcv);
  }
  QMP_sum_double_array((double *)&qcv, 2*QLA_Nc);
  push_complex_array(L, QLA_Nc, (qhmc_complex_t *)&qcv);
  return 1;
}

static int
qopqdp_squark_point(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==5);
  squark_t *q = qopqdp_squark_check(L, 1);
  int nd; get_table_len(L, 2, &nd);
  int point[nd]; get_int_array(L, 2, nd, point);
  int color = luaL_checkint(L, 3);
  double re = luaL_checknumber(L, 4);
  double im = luaL_checknumber(L, 5);
  int node = QDP_node_number(point);
  if(node==QDP_this_node) {
    int index = QDP_index(point);
    QLA_ColorVector *qcv = QDP_site_ptr_readwrite_V(q->cv, index);
    QLA_c_eq_r_plus_ir(QLA_elem_V(*qcv,color), re, im);
  }
  return 0;
}

static int
qopqdp_squark_random(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  squark_t *q = qopqdp_squark_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    sub = qopqdp_check_subset(L, 2);
  }
  QDP_V_eq_gaussian_S(q->cv, qopqdp_srs, sub);
  QLA_Real r = sqrt(0.5); // normalize to sigma^2 = 1/2
  QDP_V_eq_r_times_V(q->cv, &r, q->cv, sub);
  return 0;
}

static void
lnormalize_V(NCPROT QLA_ColorVector(*x), int i)
{
  for(int ic=0; ic<QLA_Nc; ic++) {
    QLA_Complex z = QLA_elem_V(*x,ic);
    QLA_Real n = QLA_norm2_c(z);
    if(n!=0) n = 1/sqrt(n);
    QLA_c_eq_r_times_c(QLA_elem_V(*x,ic), n, z);
  }
}

static int
qopqdp_squark_randomU1(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  squark_t *q = qopqdp_squark_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    sub = qopqdp_check_subset(L, 2);
  }
  QDP_V_eq_gaussian_S(q->cv, qopqdp_srs, sub);
  QDP_V_eq_funci(q->cv, lnormalize_V, sub);
  return 0;
}

static int
qopqdp_squark_set(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==2 || narg==3);
  squark_t *q1 = qopqdp_squark_check(L, 1);
  squark_t *q2 = qopqdp_squark_check(L, 2);
  QDP_Subset sub = QDP_all;
  if(narg!=2) {
    sub = qopqdp_check_subset(L, 3);
  }
  QDP_V_eq_V(q1->cv, q2->cv, sub);
  return 0;
}

static int
qopqdp_squark_norm2(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  squark_t *q = qopqdp_squark_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    const char *s = luaL_checkstring(L, 2);
    if(!strcmp(s,"timeslices")) sub = NULL;
    else sub = qopqdp_check_subset(L, 2);
  }
  if(sub) {
    QLA_Real nrm2;
    QDP_r_eq_norm2_V(&nrm2, q->cv, sub);
    lua_pushnumber(L, nrm2);
  } else { // timeslices
    int nt = QDP_coord_size(3);
    QLA_Real nrm2[nt];
    QDP_Subset *ts = qhmcqdp_get_timeslices();
    QDP_r_eq_norm2_V_multi(nrm2, q->cv, ts, nt);
    push_double_array(L, nt, nrm2);
  }
  return 1;
}

static int
qopqdp_squark_redot(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==2 || narg==3);
  squark_t *q1 = qopqdp_squark_check(L, 1);
  squark_t *q2 = qopqdp_squark_check(L, 2);
  QDP_Subset sub = QDP_all;
  if(narg!=2) {
    const char *s = luaL_checkstring(L, 3);
    if(!strcmp(s,"timeslices")) sub = NULL;
    else sub = qopqdp_check_subset(L, 3);
  }
  if(sub) {
    QLA_Real redot;
    QDP_r_eq_re_V_dot_V(&redot, q1->cv, q2->cv, sub);
    lua_pushnumber(L, redot);
  } else {
    int nt = QDP_coord_size(3);
    QLA_Real redot[nt];
    QDP_Subset *ts = qhmcqdp_get_timeslices();
    QDP_r_eq_re_V_dot_V_multi(redot, q1->cv, q2->cv, ts, nt);
    push_double_array(L, nt, redot);
  }
  return 1;
}

// 1: squark
// 2: squark
// 3: (optional) subset
static int
qopqdp_squark_dot(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==2 || narg==3);
  squark_t *q1 = qopqdp_squark_check(L, 1);
  squark_t *q2 = qopqdp_squark_check(L, 2);
  QDP_Subset sub = QDP_all;
  if(narg!=2) {
    const char *s = luaL_checkstring(L, 3);
    if(!strcmp(s,"timeslices")) sub = NULL;
    else sub = qopqdp_check_subset(L, 3);
  }
  if(sub) {
    QLA_Complex dot;
    QDP_c_eq_V_dot_V(&dot, q1->cv, q2->cv, sub);
    qhmc_complex_create(L, QLA_real(dot), QLA_imag(dot));
  } else { // timeslices
    int nt = QDP_coord_size(3);
    QLA_Complex dot[nt];
    QDP_Subset *ts = qhmcqdp_get_timeslices();
    QDP_c_eq_V_dot_V_multi(dot, q1->cv, q2->cv, ts, nt);
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
qopqdp_squark_combine(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==3 || narg==4);
  squark_t *qd = qopqdp_squark_check(L, 1);
  int nqs; get_table_len(L, 2, &nqs);
  squark_t *qs[nqs]; qopqdp_squark_array_check(L, 2, nqs, qs);
  int nc; get_table_len(L, 3, &nc);
  qassert(nqs==nc);
  QDP_Subset sub = QDP_all;
  if(narg>3) {
    sub = qopqdp_check_subset(L, 4);
  }
  for(int i=0; i<nqs; i++) {
    lua_pushinteger(L, i+1);
    lua_gettable(L, 3);
    if(lua_type(L,-1)==LUA_TNUMBER) {
      QLA_Real r = lua_tonumber(L, -1);
      if(i==0) QDP_V_eq_r_times_V(qd->cv, &r, qs[0]->cv, sub);
      else QDP_V_peq_r_times_V(qd->cv, &r, qs[i]->cv, sub);
    } else { // complex
      qhmc_complex_t *c = qhmc_complex_check(L, -1);
      QLA_Complex z;
      QLA_c_eq_r_plus_ir(z, c->r, c->i);
      if(i==0) QDP_V_eq_c_times_V(qd->cv, &z, qs[0]->cv, sub);
      else QDP_V_peq_c_times_V(qd->cv, &z, qs[i]->cv, sub);
    }
  }
  return 0;
}

static int
qopqdp_squark_symshift(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==4);
  squark_t *qd = qopqdp_squark_check(L, 1);
  squark_t *qs = qopqdp_squark_check(L, 2);
  gauge_t *g = qopqdp_gauge_check(L, 3);
  int mu = luaL_checkint(L, 4) - 1;
  QDP_ColorVector *tf = QDP_create_V();
  QDP_ColorVector *tb1 = QDP_create_V();
  QDP_ColorVector *tb2 = QDP_create_V();
  QDP_V_eq_sV(tf, qs->cv, QDP_neighbor[mu], QDP_forward, QDP_all);
  QDP_V_eq_Ma_times_V(tb1, g->links[mu], qs->cv, QDP_all);
  QDP_V_eq_sV(tb2, tb1, QDP_neighbor[mu], QDP_backward, QDP_all);
  QDP_V_eq_M_times_V(qd->cv, g->links[mu], tf, QDP_all);
  QDP_V_peq_V(qd->cv, tb2, QDP_all);
  QDP_destroy_V(tf);
  QDP_destroy_V(tb1);
  QDP_destroy_V(tb2);
  return 0;
}

/* ESW additions for staggered quark rephasing */

// Define a global parameter to store a phase.
// This being global is sadly unavoidable.
// It's how the gauge field rephasing is also done,
// so it's presumably a safe method to use.
int phase_squark_esw;

// Define a global parameter to store a relative location.
int *relative_loc_esw;

// Here's a function used to modify the source field with vectors.
// Based on void point_V in qdp-1.9.1/examples/sample_ks.c
void private_squark_rephase(QLA_ColorVector *s, int coords[])
{
  // Build up the phase.
  int phase = 0;
  int bits = phase_squark_esw;
  for (int i = 0; bits; i++)
    {
      // if the coordinate minus the relative coordinate is odd,
      // and if the phase flag is set, add 1.
      // if the coordinate minus the relative coordinate is odd,
      // and if the phase flag is not set, add 0.
      // if the coordinate minus the relative coordinate is even,
      // always add 0.

      if(bits&1) {
	phase += coords[i] - relative_loc_esw[i];
      }
      bits >>= 1;
    }

  if ((phase%2) == 1) // if the overall phase exponent is odd
    {
      QLA_V_eqm_V(s, s);
    }
}

static int
qopqdp_squark_rephase(lua_State *L)
{
  int narg = lua_gettop(L); // Get the number of arguments.
  qassert(narg==3);
  squark_t *qd = qopqdp_squark_check(L, 1); // Check for a staggered quark.
  int phase_flag = luaL_checkint(L, 2); // Check for the phase flag.
  int nd; get_table_len(L, 3, &nd); // Get the length of the relative coordinate.
  int r0[nd];
  get_int_array(L, 3, nd, r0);
  relative_loc_esw = r0;

  // The phase flag is defined bit-wise.
  // 0 = no phase, 1 = phase
  // (x4_sign -- time) (x3_sign) (x2_sign) (x1_sign)
  phase_squark_esw = phase_flag;

  // Apply the phase.
  QDP_V_eq_func(qd->cv, private_squark_rephase, QDP_all);

  return 0;
}

static struct luaL_Reg squark_reg[] = {
  { "__gc",     qopqdp_squark_gc },
  { "clone",    qopqdp_squark_clone },
  { "zero",     qopqdp_squark_zero },
  { "getSite",  qopqdp_squark_getSite },
  { "point",    qopqdp_squark_point },
  { "random",   qopqdp_squark_random },
  { "randomU1", qopqdp_squark_randomU1 },
  { "set",      qopqdp_squark_set },
  { "norm2",    qopqdp_squark_norm2 },
  { "Re_dot",   qopqdp_squark_redot },
  { "dot",      qopqdp_squark_dot },
  { "combine",  qopqdp_squark_combine },
  { "symshift", qopqdp_squark_symshift },
  { "rephase",  qopqdp_squark_rephase }, // ESW addition 8/23/2013
  { NULL, NULL}
};

squark_t *
qopqdp_squark_create_unset(lua_State *L)
{
  squark_t *q = lua_newuserdata(L, sizeof(squark_t));
  q->cv = QDP_create_V();
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, squark_reg);
  }
  lua_setmetatable(L, -2);
  return q;
}

squark_t *
qopqdp_squark_create(lua_State *L)
{
  squark_t *q = qopqdp_squark_create_unset(L);
  QDP_V_eq_zero(q->cv, QDP_all);
  return q;
}
