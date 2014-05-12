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
  BEGIN_ARGS;
  GET_SQUARK(q2);
  OPT_SUBSET(sub, q2->lat, QDP_all_L(q2->qlat));
  END_ARGS;
  squark_t *q1 = qopqdp_squark_create_unset(L, q2->nc, q2->lat);
  QDP_V_eq_V(q1->cv, q2->cv, sub);
  return 1;
}

static int
qopqdp_squark_zero(lua_State *L)
{
  BEGIN_ARGS;
  GET_SQUARK(q);
  OPT_SUBSET(sub, q->lat, QDP_all_L(q->qlat));
  END_ARGS;
  QDP_V_eq_zero(q->cv, sub);
  return 0;
}

// color starts at 1
static int
qopqdp_squark_point(lua_State *L)
{
  BEGIN_ARGS;
  GET_SQUARK(q);
  GET_TABLE_LEN_INDEX(nd,ip);
  GET_INT(color);
  GET_AS_COMPLEX(z);
  END_ARGS;
  int point[nd]; qhmc_get_int_array(L, ip, nd, point);
  int node = QDP_node_number_L(q->qlat, point);
  if(node==QDP_this_node) {
    int index = QDP_index_L(q->qlat, point);
    QLA_ColorVector *qcv = QDP_site_ptr_readwrite_V(q->cv, index);
    QLA_c_eq_r_plus_ir(QLA_elem_V(*qcv,color-1), z.r, z.i);
  }
  return 0;
}

static int
qopqdp_squark_random(lua_State *L)
{
  BEGIN_ARGS;
  GET_SQUARK(q);
  OPT_SUBSET(sub, q->lat, QDP_all_L(q->qlat));
  END_ARGS;
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
  BEGIN_ARGS;
  GET_SQUARK(q);
  OPT_SUBSET(sub, q->lat, QDP_all_L(q->qlat));
  END_ARGS;
  QDP_V_eq_gaussian_S(q->cv, qopqdp_srs, sub);
  QDP_V_eq_funci(q->cv, lnormalize_V, sub);
  return 0;
}

static int
qopqdp_squark_set(lua_State *L)
{
  BEGIN_ARGS;
  GET_SQUARK(q1);
  GET_SQUARK(q2);
  OPT_SUBSET(sub, q1->lat, QDP_all_L(q1->qlat));
  END_ARGS;
  QDP_V_eq_V(q1->cv, q2->cv, sub);
  return 0;
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
    sub = qopqdp_check_qsubset(L, 4, qd->lat);
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
qopqdp_squark_getSite(lua_State *L)
{
#define NC QDP_get_nc(q->cv)
  BEGIN_ARGS;
  GET_SQUARK(q);
  GET_TABLE_LEN_INDEX(nd,ip);
  END_ARGS;
  int site[nd]; qhmc_get_int_array(L, ip, nd, site);
  QLA_ColorVector qcv;
  int node = QDP_node_number_L(q->qlat, site);
  if(node==QDP_this_node) {
    int index = QDP_index_L(q->qlat, site);
    QLA_ColorVector *qcvp = QDP_site_ptr_readonly_V(q->cv, index);
    QLA_V_eq_V(&qcv, qcvp);
  } else {
    QLA_V_eq_zero(&qcv);
  }
  QMP_sum_double_array((double *)&qcv, 2*QLA_Nc);
  qhmc_push_complex_array(L, QLA_Nc, (qhmc_complex_t *)&qcv);
  return 1;
#undef NC
}

static int
qopqdp_squark_norm2(lua_State *L)
{
  BEGIN_ARGS;
  GET_SQUARK(q);
  OPT_SUBSETS(subs, ns, q->lat, QDP_all_and_empty_L(q->qlat), -1);
  END_ARGS;
  if(ns==-1) {
    QLA_Real nrm2;
    QDP_r_eq_norm2_V(&nrm2, q->cv, subs[0]);
    lua_pushnumber(L, nrm2);
  } else {
    QLA_Real nrm2[ns];
    QDP_r_eq_norm2_V_multi(nrm2, q->cv, subs, ns);
    qhmc_push_double_array(L, ns, nrm2);
  }
  return 1;
}

static int
qopqdp_squark_redot(lua_State *L)
{
  BEGIN_ARGS;
  GET_SQUARK(q1);
  GET_SQUARK(q2);
  OPT_SUBSETS(subs, ns, q1->lat, QDP_all_and_empty_L(q1->qlat), -1);
  END_ARGS;
  if(ns==-1) {
    QLA_Real redot;
    QDP_r_eq_re_V_dot_V(&redot, q1->cv, q2->cv, subs[0]);
    lua_pushnumber(L, redot);
  } else {
    QLA_Real redot[ns];
    QDP_r_eq_re_V_dot_V_multi(redot, q1->cv, q2->cv, subs, ns);
    qhmc_push_double_array(L, ns, redot);
  }
  return 1;
}

// 1: squark
// 2: squark
// 3: (optional) subset
static int
qopqdp_squark_dot(lua_State *L)
{
  BEGIN_ARGS;
  GET_SQUARK(q1);
  GET_SQUARK(q2);
  OPT_SUBSETS(subs, ns, q1->lat, QDP_all_and_empty_L(q1->qlat), -1);
  END_ARGS;
  if(ns==-1) {
    QLA_Complex dot;
    QDP_c_eq_V_dot_V(&dot, q1->cv, q2->cv, subs[0]);
    qhmc_complex_create(L, QLA_real(dot), QLA_imag(dot));
  } else {
    QLA_Complex dot[ns];
    QDP_c_eq_V_dot_V_multi(dot, q1->cv, q2->cv, subs, ns);
    qhmc_complex_t qdot[ns];
    for(int i=0; i<ns; i++) {
      qdot[i].r = QLA_real(dot[i]);
      qdot[i].i = QLA_imag(dot[i]);
    }
    qhmc_push_complex_array(L, ns, qdot);
  }
  return 1;
}

struct squark_wall_info
{
  QLA_ColorVector *value;
  int timeslice;
  int td;
};

// Value must be a pointer to a QLA_ColorVector
void
private_squark_wall(NCPROT QLA_ColorVector(*s), int coords[], void *value)
{
  struct squark_wall_info data = *(struct squark_wall_info*)value;
  if (coords[data.td] == data.timeslice) {
    QLA_V_eq_V(s, data.value);
  }
}

// Constructs a wall source. This matches the convention in MILC.
// Arguments:
// 1 - the staggered quark, as always.
// 2 - Time slice.
// 3 - An integer. 0 = even source, 1 = odd source.
// THIS IS DIFFERENT FROM THE DEFINITIONS IN GUPTA ET AL (1991).
// This is even/odd in terms of lattice sites. 
// 4 - Color
// 5 - value (real or complex)
static int
qopqdp_squark_wall(lua_State *L)
{
#define NC QDP_get_nc(q->cv)
  BEGIN_ARGS;
  GET_SQUARK(q);
  GET_INT(timeslice);
  GET_INT(evenodd);
  GET_INT(color);
  GET_AS_COMPLEX(z);
  END_ARGS;
  struct squark_wall_info data;
  QLA_ColorVector(v);
  data.value = &v;
  QLA_V_eq_zero(&v);
  QLA_c_eq_r_plus_ir(QLA_elem_V(v,color-1), z.r, z.i);
  data.timeslice = timeslice;
  data.td = q->lat->nd-1;
  QDP_V_eq_funca(q->cv, private_squark_wall, (void*)(&data), (evenodd == 0 ? QDP_even : QDP_odd));
  return 0;
#undef NC
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
#define NC nc
static void
private_squark_rephase(NCPROT QLA_ColorVector *s, int coords[])
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
#undef NC

static int
qopqdp_squark_rephase(lua_State *L)
{
  int narg = lua_gettop(L); // Get the number of arguments.
  qassert(narg==3);
  squark_t *qd = qopqdp_squark_check(L, 1); // Check for a staggered quark.
  int phase_flag = luaL_checkint(L, 2); // Check for the phase flag.
  int nd; get_table_len(L, 3, &nd); // Get the length of the relative coordinate.
  int r0[nd];
  qhmc_get_int_array(L, 3, nd, r0);
  relative_loc_esw = r0;

  // The phase flag is defined bit-wise.
  // 0 = no phase, 1 = phase
  // (x4_sign -- time) (x3_sign) (x2_sign) (x1_sign)
  phase_squark_esw = phase_flag;

  // Apply the phase.
  QDP_V_eq_func(qd->cv, private_squark_rephase, QDP_all);

  return 0;
}

static int
qopqdp_squark_symshift(lua_State *L)
{
#define NC QDP_get_nc(qs->cv)
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
#undef NC
}

// ESW addition for computing baryon correlators.
// 1: first quark
// 2: array of Nc-1 quarks
// 3: (optional) sublattice
static int
qopqdp_squark_epscontr(lua_State *L)
{
#define NC nc
  BEGIN_ARGS;
  GET_SQUARK(q1);
  GET_TABLE_LEN_INDEX(nc1,iq);
  int nc = nc1 + 1;
  squark_t *qs[nc]; qs[0]=q1; qopqdp_squark_array_check(L, iq, nc1, qs+1);
  OPT_SUBSETS(subs, ns, qs[0]->lat, QDP_all_and_empty_L(qs[0]->qlat), -1);
  END_ARGS;
  QDP_Subset sub = subs[0];
  if(ns>1) sub = QDP_all_L(qs[0]->qlat);
  // Create a color matrix, insert fermions into it.
  QDP_ColorMatrix *mat = QDP_create_M_L(qs[0]->qlat);
  QDP_Complex *determ = QDP_create_C_L(qs[0]->qlat);
  for(int i=0; i<nc; i++) {
    QDP_M_eq_colorvec_V(mat, qs[i]->cv, i, sub);
  }
  QDP_C_eq_det_M(determ, mat, sub);
  if(ns==-1) {
    QLA_Complex val;
    QDP_c_eq_sum_C(&val, determ, sub);
    qhmc_complex_create(L, QLA_real(val), QLA_imag(val));
  } else {
    QLA_Complex val[ns];
    QDP_c_eq_sum_C_multi(val, determ, subs, ns);
    qhmc_complex_t qval[ns];
    for(int i=0; i<ns; i++) {
      qval[i].r = QLA_real(val[i]);
      qval[i].i = QLA_imag(val[i]);
    }
    qhmc_push_complex_array(L, ns, qval);
  }
  QDP_destroy_M(mat);
  QDP_destroy_C(determ);
  return 1;
#undef NC
}

static int
qopqdp_squark_s4_ferm_observables(lua_State *L)
{
#define NC QDP_get_nc(q->cv)
  int narg = lua_gettop(L);
  qassert(narg==3);
  squark_t *q = qopqdp_squark_check(L, 1);
  squark_t *q2 = qopqdp_squark_check(L, 2);
  gauge_t *g = qopqdp_gauge_check(L, 3);
  int i;
  int r0[4]; // for phases.
  
  QDP_Subset* dirsubset;
  
  QLA_Complex interm[2];
  qhmc_complex_t pbp_e[6], pbp_o[6];
  // 0 - pbp on even and odd sites.
  // 1 - add t one-link pbp if x even/odd
  // 2 - add t one-link pbp if t even/odd
  // 3 - add x one-link pbp if x even/odd
  // 4 - add y one-link pbp if y even/odd
  // 5 - add z one-link pbp if z even/odd
  
  // Because we need it.
  QDP_ColorVector* temp_vec;
  QDP_Complex* temp_cmplx;
  
  // Zero out arrays.
  for (i = 0; i < 6; i++)
  {
    pbp_e[i].r = pbp_e[i].i = pbp_o[i].r = pbp_o[i].i = 0;
  }
  
  // Get even, odd pbp.
  QDP_c_eq_V_dot_V(&interm[0], q->cv, q->cv, QDP_even);
  QDP_c_eq_V_dot_V(&interm[1], q->cv, q->cv, QDP_odd);
  pbp_e[0].r = QLA_real(interm[0]);
  pbp_e[0].i = QLA_imag(interm[0]);
  pbp_o[0].r = QLA_real(interm[1]);
  pbp_o[0].i = QLA_imag(interm[1]);
  
  // Generate one-link separated.
  temp_vec = QDP_create_V();
  temp_cmplx = QDP_create_C();
  
  // Time one link.
  QDP_V_eq_M_times_sV(temp_vec, g->links[3], q->cv, QDP_neighbor[3], QDP_forward, QDP_all);
  QDP_C_eq_V_dot_V(temp_cmplx, q2->cv, temp_vec, QDP_all);
  
  dirsubset = qhmcqdp_get_eodir(q->lat, 0); // x even/odd
  QDP_c_eq_sum_C_multi(interm, temp_cmplx, dirsubset, 2);
  pbp_e[1].r = QLA_real(interm[0]);
  pbp_e[1].i = QLA_imag(interm[0]);
  pbp_o[1].r = QLA_real(interm[1]);
  pbp_o[1].i = QLA_imag(interm[1]);
  
  dirsubset = qhmcqdp_get_eodir(q->lat, 3); // t even/odd
  QDP_c_eq_sum_C_multi(interm, temp_cmplx, dirsubset, 2);
  pbp_e[2].r = QLA_real(interm[0]);
  pbp_e[2].i = QLA_imag(interm[0]);
  pbp_o[2].r = QLA_real(interm[1]);
  pbp_o[2].i = QLA_imag(interm[1]);
  
  // Prepare for phases.
  r0[0] = r0[1] = r0[2] = r0[3] = 0;
  relative_loc_esw = r0;
  
  // X one link. Don't forget phase factors!
  QDP_V_eq_M_times_sV(temp_vec, g->links[0], q->cv, QDP_neighbor[0], QDP_forward, QDP_all);
  phase_squark_esw = 0x8; // time
  QDP_V_eq_func(temp_vec, private_squark_rephase, QDP_all);
  QDP_C_eq_V_dot_V(temp_cmplx, q2->cv, temp_vec, QDP_all);
  dirsubset = qhmcqdp_get_eodir(q->lat, 0); // x even/odd
  QDP_c_eq_sum_C_multi(interm, temp_cmplx, dirsubset, 2);
  pbp_e[3].r = QLA_real(interm[0]);
  pbp_e[3].i = QLA_imag(interm[0]);
  pbp_o[3].r = QLA_real(interm[1]);
  pbp_o[3].i = QLA_imag(interm[1]);
  
  // Y one link.
  QDP_V_eq_M_times_sV(temp_vec, g->links[1], q->cv, QDP_neighbor[1], QDP_forward, QDP_all);
  phase_squark_esw = 0x9; // time, x
  QDP_V_eq_func(temp_vec, private_squark_rephase, QDP_all);
  QDP_C_eq_V_dot_V(temp_cmplx, q2->cv, temp_vec, QDP_all);
  dirsubset = qhmcqdp_get_eodir(q->lat, 1); // y even/odd
  QDP_c_eq_sum_C_multi(interm, temp_cmplx, dirsubset, 2);
  pbp_e[4].r = QLA_real(interm[0]);
  pbp_e[4].i = QLA_imag(interm[0]);
  pbp_o[4].r = QLA_real(interm[1]);
  pbp_o[4].i = QLA_imag(interm[1]);
  
  // Z one link.
  QDP_V_eq_M_times_sV(temp_vec, g->links[2], q->cv, QDP_neighbor[2], QDP_forward, QDP_all);
  phase_squark_esw = 0xB; // time, x, y
  QDP_V_eq_func(temp_vec, private_squark_rephase, QDP_all);
  QDP_C_eq_V_dot_V(temp_cmplx, q2->cv, temp_vec, QDP_all);
  dirsubset = qhmcqdp_get_eodir(q->lat, 2); // z even/odd
  QDP_c_eq_sum_C_multi(interm, temp_cmplx, dirsubset, 2);
  pbp_e[5].r = QLA_real(interm[0]);
  pbp_e[5].i = QLA_imag(interm[0]);
  pbp_o[5].r = QLA_real(interm[1]);
  pbp_o[5].i = QLA_imag(interm[1]);
  
  QDP_destroy_V(temp_vec);
  QDP_destroy_C(temp_cmplx);
  
  qhmc_push_complex_array(L, 6, pbp_e);
  qhmc_push_complex_array(L, 6, pbp_o);
  
  return 2;
#undef NC
}

static struct luaL_Reg squark_reg[] = {
  { "__gc",        qopqdp_squark_gc },
  { "clone",       qopqdp_squark_clone },
  { "zero",        qopqdp_squark_zero },
  { "point",       qopqdp_squark_point },
  { "random",      qopqdp_squark_random },
  { "randomU1",    qopqdp_squark_randomU1 },
  { "set",         qopqdp_squark_set },
  { "combine",     qopqdp_squark_combine },
  { "getSite",     qopqdp_squark_getSite },
  { "norm2",       qopqdp_squark_norm2 },
  { "reDot",       qopqdp_squark_redot },
  { "Re_dot",      qopqdp_squark_redot },
  { "dot",         qopqdp_squark_dot },
  { "wall",        qopqdp_squark_wall }, // ESW addition 10/8/2013
  { "rephase",     qopqdp_squark_rephase }, // ESW addition 8/23/2013
  { "symshift",    qopqdp_squark_symshift },
  { "epsContract", qopqdp_squark_epscontr }, // ESW addition 11/15/2013
  { "s4Ferm",      qopqdp_squark_s4_ferm_observables }, // ESW addition 12/18/2013
  { NULL, NULL}
};

squark_t *
qopqdp_squark_create_unset(lua_State *L, int nc, lattice_t *lat)
{
#define NC nc
  if(nc==0) nc = QOPQDP_DEFAULTNC;
  if(lat==NULL) lat = qopqdp_get_default_lattice(L);
  squark_t *q = lua_newuserdata(L, sizeof(squark_t));
  q->cv = QDP_create_V_L(lat->qlat);
  q->lat = lat;
  q->qlat = lat->qlat;
  q->nc = nc;
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, squark_reg);
  }
  lua_setmetatable(L, -2);
  return q;
}

squark_t *
qopqdp_squark_create(lua_State *L, int nc, lattice_t *lat)
{
  squark_t *q = qopqdp_squark_create_unset(L, nc, lat);
  QDP_V_eq_zero(q->cv, QDP_all);
  return q;
}
