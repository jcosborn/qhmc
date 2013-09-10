#include <string.h>
#include <math.h>
#include "qhmc_qopqdp_common.h"

static char *fmtname = "qopqdp.force";

force_t *
qopqdp_force_check(lua_State *L, int idx)
{
  luaL_checkudata(L, idx, fmtname);
  force_t *g = lua_touserdata(L, idx);
#if 0
  int hasmt = lua_getmetatable(L, idx);
  qassert(hasmt==1);
  luaL_getmetatable(L, fmtname);
  int eq = lua_equal(L, -1, -2);
  qassert(eq==1);
  lua_pop(L, 2);
#endif
  return g;
}

void
qopqdp_force_array_check(lua_State *L, int idx, int n, force_t *g[n])
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qassert(lua_objlen(L, idx)==n);
  lua_pushvalue(L, idx); // make copy for indexing convenience
  for(int i=0; i<n; i++) {
    lua_pushnumber(L, i+1);
    lua_gettable(L, -2);
    g[i] = qopqdp_force_check(L, -1);
    lua_pop(L, 1);
  }
  lua_pop(L, 1);
}

static void
qopqdp_force_free(lua_State *L, int idx)
{
  force_t *g = qopqdp_force_check(L, idx);
  for(int i=0; i<g->nd; i++) {
    QDP_destroy_M(g->force[i]);
  }
}

static int
qopqdp_force_gc(lua_State *L)
{
  qopqdp_force_free(L, -1);
  return 0;
}

static int
qopqdp_force_clone(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==1 || narg==2);
  force_t *f1 = qopqdp_force_create(L);
  force_t *f2 = qopqdp_force_check(L, 1);
  QDP_Subset sub = QDP_all;
  if(narg!=1) {
    sub = qopqdp_check_subset(L, 2);
  }
  for(int i=0; i<f1->nd; i++) {
    QDP_M_eq_M(f1->force[i], f2->force[i], sub);
  }
  return 1;
}

static int
qopqdp_force_zero(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  force_t *f = qopqdp_force_check(L, -1);
  for(int i=0; i<f->nd; i++) {
    QDP_M_eq_zero(f->force[i], QDP_all);
  }
  return 0;
}

static int
qopqdp_force_set(lua_State *L)
{
  qassert(lua_gettop(L)==2);
  force_t *f1 = qopqdp_force_check(L, 1);
  force_t *f2 = qopqdp_force_check(L, 2);
  for(int i=0; i<f1->nd; i++) {
    QDP_M_eq_M(f1->force[i], f2->force[i], QDP_all);
  }
  return 0;
}

static void
randforce(NCPROT QLA_ColorMatrix(*m), int i, void *args)
{
  QLA_RandomState *s = (QLA_RandomState*)args + i;
  QLA_Real s2 = 0.70710678118654752440;  // sqrt(1/2)
  QLA_Real s3 = 0.57735026918962576450;  // sqrt(1/3)
  QLA_Real r3 = s2*QLA_gaussian(s);
  QLA_Real r8 = s2*s3*QLA_gaussian(s);
  QLA_c_eq_r_plus_ir(QLA_elem_M(*m,0,0), 0, r3+r8);
  QLA_c_eq_r_plus_ir(QLA_elem_M(*m,1,1), 0, -r3+r8);
  QLA_c_eq_r_plus_ir(QLA_elem_M(*m,2,2), 0, -2*r8);
  QLA_Real r01 = s2*QLA_gaussian(s);
  QLA_Real r02 = s2*QLA_gaussian(s);
  QLA_Real r12 = s2*QLA_gaussian(s);
  QLA_Real i01 = s2*QLA_gaussian(s);
  QLA_Real i02 = s2*QLA_gaussian(s);
  QLA_Real i12 = s2*QLA_gaussian(s);
  QLA_c_eq_r_plus_ir(QLA_elem_M(*m,0,1),  r01, i01);
  QLA_c_eq_r_plus_ir(QLA_elem_M(*m,1,0), -r01, i01);
  QLA_c_eq_r_plus_ir(QLA_elem_M(*m,0,2),  r02, i02);
  QLA_c_eq_r_plus_ir(QLA_elem_M(*m,2,0), -r02, i02);
  QLA_c_eq_r_plus_ir(QLA_elem_M(*m,1,2),  r12, i12);
  QLA_c_eq_r_plus_ir(QLA_elem_M(*m,2,1), -r12, i12);
}

static int
qopqdp_force_random(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  force_t *f = qopqdp_force_check(L, 1);
  if(QLA_Nc==3) { // use MILC conventions
    QLA_RandomState *s = QDP_expose_S(qopqdp_srs);
    for(int i=0; i<f->nd; i++) {
      QDP_M_eq_funcia(f->force[i], randforce, s, QDP_all);
      //{ QLA_Real s=1.1; QDP_M_eq_r_times_M(f->force[i], &s, f->force[i], QDP_all); }
    }
    QDP_reset_S(qopqdp_srs);
  } else {
    QDP_ColorMatrix *m = QDP_create_M();
    for(int i=0; i<f->nd; i++) {
      QDP_M_eq_gaussian_S(m, qopqdp_srs, QDP_all);
      QDP_M_eq_antiherm_M(f->force[i], m, QDP_all);
    }
    QDP_destroy_M(m);
  }
  return 0;
}

// 1: force
// 2: number or table of numbers
static int
qopqdp_force_scale(lua_State *L)
{
  qassert(lua_gettop(L)==2);
  force_t *f = qopqdp_force_check(L, 1);
  int nd = f->nd;
  double sa[nd];
  if(lua_type(L,2)==LUA_TTABLE) {
    get_double_array(L, 2, nd, sa);
  } else {
    double s = luaL_checknumber(L, 2);
    for(int i=0; i<nd; i++) sa[i] = s;
  }
  for(int i=0; i<nd; i++) {
    QLA_Real s = sa[i];
    if(s!=1) QDP_M_eq_r_times_M(f->force[i], &s, f->force[i], QDP_all);
  }
  return 0;
}

static int
qopqdp_force_norm2(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  force_t *f = qopqdp_force_check(L, 1);
  QLA_Real nrm = 0;
  double tn[f->nd];
  for(int i=0; i<f->nd; i++) {
    QLA_Real t;
    QDP_r_eq_norm2_M(&t, f->force[i], QDP_all);
    nrm += t;
    tn[i] = t;
  }
  lua_pushnumber(L, nrm);
  push_double_array(L, f->nd, tn);
  return 2;
}

static int
qopqdp_force_infnorm(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  force_t *f = qopqdp_force_check(L, 1);
  QLA_Real nrm = 0;
  for(int i=0; i<f->nd; i++) {
    QLA_Real t = infnorm_M(f->force[i], QDP_all);
    if(t>nrm) nrm = t;
  }
  lua_pushnumber(L, nrm);
  return 1;
}

static int
qopqdp_force_derivForce(lua_State *L)
{
  qassert(lua_gettop(L)==2);
  force_t *f = qopqdp_force_check(L, 1);
  gauge_t *g = qopqdp_gauge_check(L, 2);
  QDP_ColorMatrix *t = QDP_create_M();
  for(int mu=0; mu<4; mu++) {
    QDP_M_eq_M_times_Ma(t, g->links[mu], f->force[mu], QDP_all);
    // factor of -2 for GL -> U
    QDP_M_eq_r_times_M(t, (QLA_Real[1]){-2}, t, QDP_all);
    QDP_M_eq_antiherm_M(f->force[mu], t, QDP_all);
  }
  //info->final_flop += (4.*(198+24+18))*QDP_sites_on_node;
  QDP_destroy_M(t);
  return 0;
}

static int
qopqdp_force_update(lua_State *L)
{
  qassert(lua_gettop(L)==3);
  force_t *f1 = qopqdp_force_check(L, 1);
  force_t *f2 = qopqdp_force_check(L, 2);
  QLA_Real eps = luaL_checknumber(L, 3);
  for(int i=0; i<f1->nd; i++) {
    QDP_M_peq_r_times_M(f1->force[i], &eps, f2->force[i], QDP_all);
  }
  return 0;
}

static int
qopqdp_force_time(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  force_t *f = qopqdp_force_check(L, 1);
  lua_pushnumber(L, f->time);
  return 1;
}

static int
qopqdp_force_flops(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  force_t *f = qopqdp_force_check(L, 1);
  lua_pushnumber(L, f->flops);
  return 1;
}

#if 0
static void
get_gen(QLA_ColorMatrix *g, int ig)
{
  QLA_M_eq_zero(g);
  if(ig==0) {
    for(int i=0; i<QLA_Nc; i++) {
      QLA_c_eq_r(QLA_elem_M(*g,i,i), 1);
    }
  } else {
    int c = (int) floor(sqrt((double)ig));
    int t = ig - c*c;
    int r = t/2;
    if(2*r!=t) { int s = c; c = r; r = s; }
    if(r<c) {
      QLA_c_eq_r(QLA_elem_M(*g,r,c), 1);
      QLA_c_eq_r(QLA_elem_M(*g,c,r), 1);
    } else if(c<r) {
      QLA_c_eq_r_plus_ir(QLA_elem_M(*g,r,c), 0, 1);
      QLA_c_eq_r_plus_ir(QLA_elem_M(*g,c,r), 0, -1);
    } else {
      QLA_Real s = sqrt(2./(r*(r+1)));
      for(int i=0; i<r; i++) {
	QLA_c_eq_r(QLA_elem_M(*g,i,i), s);
      }
      QLA_c_eq_r(QLA_elem_M(*g,r,r), -s*r);
    }
  }
}

void
get_local_force(QLA_ColorMatrix *f, QLA_ColorMatrix *g0,
		void (*update_link)(QLA_ColorMatrix *g, QLA_ColorMatrix *g0,
				    QLA_ColorMatrix *dg, double eps),
		double (*act)(QLA_ColorMatrix *g, void *), void *args, double eps)
{
  QLA_M_eq_zero(f);
  double s0 = act(g0, args);
  for(int ig=1; ig<=8; ig++) {
    QLA_ColorMatrix gen, g;
    get_gen(&gen, ig);
    update_link(&g, g0, &gen, eps);
    double s = act(&g, args);
    QLA_Real d = -0.5*(s-s0)/eps;
    QLA_Complex z;
    QLA_c_eq_r_plus_ir(z, 0, d);
    QLA_M_peq_c_times_M(f, &z, &gen);
  }
}

static void
dump_M(QLA_ColorMatrix *m)
{
  for(int i=0; i<QLA_Nc; i++) {
    for(int j=0; j<QLA_Nc; j++) {
      QLA_Complex *z = &QLA_elem_M(*m, i, j);
      printf("%i %i\t%g\t%g\n", i, j, QLA_real(*z), QLA_imag(*z));
    }
  }
}

static void
compare_M(QLA_ColorMatrix *m1, QLA_ColorMatrix *m2)
{
  QLA_ColorMatrix d;
  QLA_Real n;
  QLA_M_eq_M_minus_M(&d, m1, m2);
  QLA_r_eq_norm2_M(&n, &d);
  if(n>1e-10) {
    printf("mat diff norm2 = %g\n", n);
    printf("mat1:\n");
    dump_M(m1);
    printf("mat2:\n");
    dump_M(m2);
  }
}

void
check_force(force_t *f, gauge_t *g, double (*act)(gauge_t *g, void *), void *args)
{
  QLA_Complex ieps;
  double eps = 1e-6;
  QLA_c_eq_r_plus_ir(ieps, 0, eps);
  double ao = act(g, args);
  for(int mu=0; mu<f->nd; mu++) {
    QLA_ColorMatrix *qf = QDP_expose_M(f->force[mu]);
    for(int idx=0; idx<QDP_sites_on_node; idx++) {
      QLA_ColorMatrix gs, ff;
      QLA_ColorMatrix *qg = QDP_expose_M(g->links[mu]);
      QLA_M_eq_M(&gs, qg+idx);
      QDP_reset_M(g->links[mu]);
#if 1
      QLA_ColorMatrix ls;
      QLA_ColorMatrix *ql = QDP_expose_M(g->lie[mu]);
      QLA_M_eq_M(&ls, ql+idx);
      QDP_reset_M(g->lie[mu]);
#endif
      QLA_M_eq_zero(&ff);
      for(int ig=1; ig<=8; ig++) {
	QLA_ColorMatrix gen, eg, eeg;
	get_gen(&gen, ig);
#if 0
	QLA_M_eq_c_times_M(&eg, &ieps, &gen);
	QLA_M_eq_exp_M(&eeg, &eg);
	qg = QDP_expose_M(g->links[mu]);
	QLA_M_eq_M_times_M(qg+idx, &eeg, &gs);
	QDP_reset_M(g->links[mu]);
#else
	qg = QDP_expose_M(g->links[mu]);
	ql = QDP_expose_M(g->lie[mu]);
	QLA_M_eq_M(ql+idx, &ls);
	QLA_M_peq_r_times_M(ql+idx, &eps, &gen);
	QLA_M_eq_i_M(&eg, ql+idx);
	QLA_M_eq_exp_M(qg+idx, &eg);
	QDP_reset_M(g->links[mu]);
	QDP_reset_M(g->lie[mu]);
#endif
	double an = act(g, args);
	QLA_Real f = 0.5*(an-ao)/eps;
	QLA_Complex z;
	QLA_c_eq_r_plus_ir(z, 0, f);
	QLA_M_eq_c_times_M(&eg, &z, &gen);
	QLA_M_peq_M(&ff, &eg);
      }
      qg = QDP_expose_M(g->links[mu]);
      QLA_M_eq_M(qg+idx, &gs);
      QDP_reset_M(g->links[mu]);
#if 1
      ql = QDP_expose_M(g->lie[mu]);
      QLA_M_eq_M(ql+idx, &ls);
      QDP_reset_M(g->lie[mu]);
#endif
      compare_M(qf+idx, &ff);
    }
    QDP_reset_M(f->force[mu]);
  }
}
#endif

static struct luaL_Reg force_reg[] = {
  { "__gc",       qopqdp_force_gc },
  { "clone",      qopqdp_force_clone },
  { "zero",       qopqdp_force_zero },
  { "set",        qopqdp_force_set },
  { "random",     qopqdp_force_random },
  { "scale",      qopqdp_force_scale },
  { "norm2",      qopqdp_force_norm2 },
  { "infnorm",    qopqdp_force_infnorm },
  { "derivForce", qopqdp_force_derivForce },
  { "update",     qopqdp_force_update },
  { "time",       qopqdp_force_time },
  { "flops",      qopqdp_force_flops },
  { NULL, NULL}
};

force_t *
qopqdp_force_create(lua_State* L)
{
  int nd = QDP_ndim();
  force_t *f = lua_newuserdata(L, sizeof(force_t)+nd*sizeof(QDP_ColorMatrix*));
  f->nd = nd;
  for(int i=0; i<nd; i++) {
    f->force[i] = QDP_create_M();
    QDP_M_eq_zero(f->force[i], QDP_all);
  }
  if(luaL_newmetatable(L, fmtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, force_reg);
  }
  lua_setmetatable(L, -2);
  return f;
}
