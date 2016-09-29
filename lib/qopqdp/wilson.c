#include <string.h>
#include "qhmc_qopqdp_common.h"
#ifdef HAVE_NC1
#include <qdp_d1.h>
#include <qop_d1.h>
#endif
#ifdef HAVE_NC2
#include <qdp_d2.h>
#include <qop_d2.h>
#endif
#ifdef HAVE_NC3
#include <qdp_d3.h>
#include <qop_d3.h>
#endif

static char *mtname = "qopqdp.wilson";

#define kappa(w,m) (0.5/(1.+3.*((w)->coeffs.aniso)+(m)))

wilson_t *
qopqdp_wilson_check(lua_State *L, int idx)
{
  luaL_checkudata(L, idx, mtname);
  wilson_t *w = lua_touserdata(L, idx);
  return w;
}

static void
qopqdp_wilson_free(lua_State *L, int idx)
{
  wilson_t *w = qopqdp_wilson_check(L, idx);
  if(w->fl) {
    QOP_wilson_destroy_L(w->fl);
    w->fl = NULL;
  }
  if(w->ffl) {
    QOP_F_wilson_destroy_L(w->ffl);
    w->ffl = NULL;
  }
}

static int
qopqdp_wilson_gc(lua_State *L)
{
  qopqdp_wilson_free(L, -1);
  return 0;
}

static void
qopqdp_wilson_set_opts(void)
{
  QOP_opt_t opt[3];
  opt[0].tag = "st";
  opt[0].value = 3;
  opt[1].tag = "ns";
  opt[1].value = 8;
  opt[2].tag = "nm";
  opt[2].value = 8;
  QOP_wilson_invert_set_opts(opt, 3);
  opt[0].tag = "cg";
  opt[0].value = 0;
  QOP_wilson_invert_set_opts(opt, 1);
  //opt[0].tag = "fnmat_src_min";
  //opt[0].value = 0;
  //QOP_wilson_force_set_opts(opt, 1);
#ifdef __bg__
  //QDP_set_block_size(64);
#else
  //QDP_set_block_size(1024);
#endif
}

#if 0
static void
rephase(QDP_ColorMatrix *g[], int nd)
{
  QOP_Complex bcphase[nd];
  QOP_bc_t bc;
  bc.phase = bcphase;
  for(int i=0; i<nd; i++) {
    bcphase[i].re = 1;
    bcphase[i].im = 0;
  }
  bcphase[nd-1].re = -1;
  int signmask[4] = {0,0,0,0};
  QOP_staggered_sign_t ss;
  ss.signmask = signmask;
  int r0[4] = {0,0,0,0};
  QOP_GaugeField *qg = QOP_create_G_from_qdp(g);
  QOP_rephase_G(qg, r0, &bc, &ss);
  QOP_extract_G_to_qdp(g, qg);
  QOP_destroy_G(qg);
}
#endif

#define set_phases(nd)	 \
  int r0[4] = {0,0,0,0};   \
  QOP_Complex bcphase[nd]; \
  QOP_bc_t bc; \
  bc.phase = bcphase; \
  int signmask[4] = {0,0,0,0}; \
  QOP_staggered_sign_t ss; \
  ss.signmask = signmask; \
  set_phases_func(r0, &bc, &ss, nd)

static void
set_phases_func(int *r0, QOP_bc_t *bc, QOP_staggered_sign_t *ss, int nd)
{
  for(int i=0; i<nd; i++) {
    bc->phase[i].re = 1;
    bc->phase[i].im = 0;
  }
#ifndef QHMC_REPRO_UNIFORM
  bc->phase[nd-1].re = -1;
#endif
}

#if 0
static void
rephase_G(QOP_GaugeField *g, int nd)
{
  set_phases(nd);
  QOP_rephase_G(g, r0, &bc, &ss);
}
#endif

static void
rephase_F(gauge_t *f)
{
  set_phases(f->nd);
  QOP_rephase_G_qdp(f->links, r0, &bc, &ss);
}

static void
wilson_set(wilson_t *w, int prec)
{
#define NC QDP_get_nc(g->links[0])
  if((prec==1&&w->ffl)||(prec==2&&w->fl)) {
    return;
  }
  gauge_t *g = w->g;
  int nd = g->nd;
  QOP_Complex bcphase[nd];
  QOP_bc_t bc;
  bc.phase = bcphase;
  for(int i=0; i<nd; i++) {
    bcphase[i].re = 1;
    bcphase[i].im = 0;
  }
#ifndef QHMC_REPRO_UNIFORM
  bcphase[nd-1].re = -1;
#endif
  int signmask[4] = {0,0,0,0};
  QOP_staggered_sign_t ss;
  ss.signmask = signmask;
  int r0[4] = {0,0,0,0};
  QOP_info_t info;
  if(prec==1) {
    QDP_F_ColorMatrix *flinks[nd];
    for(int i=0; i<nd; i++) {
      flinks[i] = QDP_F_create_M();
      QDP_FD_M_eq_M(flinks[i], g->links[i], QDP_all);
    }
    QOP_F_GaugeField *qg = QOP_F_create_G_from_qdp(flinks);
    QOP_F_rephase_G(qg, r0, &bc, &ss);
    w->ffl = QOP_F_wilson_create_L_from_G(&info, &w->coeffs, qg);
    QOP_F_destroy_G(qg);
    for(int i=0; i<nd; i++) {
      QDP_F_destroy_M(flinks[i]);
    }
  } else {
    QOP_GaugeField *qg = QOP_create_G_from_qdp(g->links);
    QOP_rephase_G(qg, r0, &bc, &ss);
    w->fl = QOP_wilson_create_L_from_G(&info, &w->coeffs, qg);
    QOP_destroy_G(qg);
  }
  w->time = info.final_sec;
  w->flops = info.final_flop;
#undef NC
}

// 1: wilson
// 2: gauge
// 3: (optional) table of coefficients (clov_s, clov_t, aniso(=nu/xi0))
// 3/4: (optional) precision
static int
qopqdp_wilson_set(lua_State *L)
{
  int nargs = lua_gettop(L);
  qassert(nargs>=2 && nargs<=4);
  wilson_t *w = qopqdp_wilson_check(L, 1);
  gauge_t *g = qopqdp_gauge_check(L, 2);
  double prec = 0;
  if(lua_type(L,3)==LUA_TTABLE) {
    tableLoopKeys(L, 3) {
      //printf0("%s\n", tableLoopKey);
      tableLoopKeysIfKeySetDouble(L, "clov_s", w->coeffs.clov_s);
      else tableLoopKeysIfKeySetDouble(L, "clov_t", w->coeffs.clov_t);
      else tableLoopKeysIfKeySetDouble(L, "aniso", w->coeffs.aniso);
      tableLoopKeysEnd(L);
    }
    prec = luaL_optint(L, 4, 0);
  } else {
    prec = luaL_optint(L, 3, 0);
  }
  if(w->fl) {
    QOP_wilson_destroy_L(w->fl);
    w->fl = NULL;
  }
  if(w->ffl) {
    QOP_F_wilson_destroy_L(w->ffl);
    w->ffl = NULL;
  }
  qopqdp_wilson_set_opts();
  w->g = g;
  if(prec==1) wilson_set(w, 1);
  if(prec==2) wilson_set(w, 2);
  return 0;
}

// args: wilson, table{ level=table{ key=value, ... } }
static int
qopqdp_wilson_mg_set(lua_State *L)
{
  int nargs = lua_gettop(L);
  qassert(nargs==2);
  wilson_t *w = qopqdp_wilson_check(L, 1);
  if(!w->mg) {
    w->mg = QOP_wilsonMgNew();
  }

  luaL_checktype(L, 2, LUA_TTABLE);
  for(int l=-2; l<10; l++) { // mg setup levels
    lua_pushnumber(L, l);
    lua_gettable(L, 2);
    if(!lua_isnil(L, -1)) { // have option table at this level
      lua_pushnil(L);  /* first key */
      while (lua_next(L, -2) != 0) {
	/* uses 'key' (at index -2) and 'value' (at index -1) */
	lua_pushvalue(L, -2);
	const char *key = luaL_checkstring(L, -1);
	printf0("%i  %s", l, key);
	lua_pop(L, 1);
	if(lua_istable(L, -1)) {
	  int n = lua_objlen(L, -1);
	  double val[n];
	  for(int i=0; i<n; i++) {
	    lua_pushnumber(L, i+1);
	    lua_gettable(L, -2);
	    val[i] = luaL_checknumber(L, -1);
	    printf0("  %g", val[i]);
	    lua_pop(L, 1);
	  }
	  printf0("\n");
	  QOP_wilsonMgSetArray(w->mg, l, (char*)key, val, n);
	} else {
	  double val = luaL_checknumber(L, -1);
	  printf0("  %g\n", val);
	  QOP_wilsonMgSet(w->mg, l, (char*)key, val);
	}
	/* removes 'value'; keeps 'key' for next iteration */
	lua_pop(L, 1);
      }
    }
    lua_pop(L, 1);
  }
  return 0;
}

// args: wilson
static int
qopqdp_wilson_mg_setup(lua_State *L)
{
  int nargs = lua_gettop(L);
  qassert(nargs==1);
  wilson_t *w = qopqdp_wilson_check(L, 1);
  wilson_set(w, 1);
  QOP_wilsonMgSet(w->mg, -1, "nc", w->g->nc);
#ifdef QOP_wilsonMgSetLinks
  QOP_F_wilsonMgSetLinks(w->mg, w->ffl);
#else
  QOP_wilsonMgSetLinks(w->mg, w->ffl);
#endif
  QOP_wilsonMgSetup(w->mg);
  return 0;
}

static int
qopqdp_wilson_D(lua_State *L)
{
  BEGIN_ARGS;
  GET_WILSON(w);
  GET_QOPQDP_DFERMION(qd);
  GET_QOPQDP_DFERMION(qs);
  GET_DOUBLE(mass);
  OPT_EVENODD(eod, QOP_EVENODD);
  OPT_EVENODD(eos, QOP_EVENODD);
  END_ARGS;
  wilson_set(w, 2);
  QOP_wilson_dslash_qdp(NULL, w->fl, kappa(w,mass), 1,
			qd->field, qs->field, eod, eos);
  return 0;
}

static int
qopqdp_wilson_Ddag(lua_State *L)
{
  BEGIN_ARGS;
  GET_WILSON(w);
  GET_QOPQDP_DFERMION(qd);
  GET_QOPQDP_DFERMION(qs);
  GET_DOUBLE(mass);
  OPT_EVENODD(eod, QOP_EVENODD);
  OPT_EVENODD(eos, QOP_EVENODD);
  END_ARGS;
  wilson_set(w, 2);
  QOP_wilson_dslash_qdp(NULL, w->fl, kappa(w,mass), -1,
			qd->field, qs->field, eod, eos);
  return 0;
}

static int
qopqdp_wilson_precD(lua_State *L)
{
#define NC QDP_get_nc(qd->field)
  BEGIN_ARGS;
  GET_WILSON(w);
  GET_QOPQDP_DFERMION(qd);
  GET_QOPQDP_DFERMION(qs);
  GET_DOUBLE(mass);
  END_ARGS;
  wilson_set(w, 2);
  double k = kappa(w,mass);
#if 0
  QOP_wilson_dslash_qdp(NULL, w->fl, k, 1, qd->df, qs->df, QOP_ODD, QOP_EVEN);
  QOP_wilson_dslash_qdp(NULL, w->fl, k, 1, qd->df, qd->df, QOP_EVEN, QOP_ODD);
  QLA_Real s = -4*k*k;
  QDP_D_eq_r_times_D_plus_D(qd->df, &s, qd->df, qs->df, QDP_even);
#else
  QDP_Lattice *lat = qs->qlat;
  QDP_DiracFermion *t1 = QDP_create_D_L(lat);
  QDP_DiracFermion *t2 = QDP_create_D_L(lat);
  QOP_wilson_diaginv_qdp(NULL, w->fl, k, t1, qs->field, QOP_EVEN);
  QOP_wilson_dslash_qdp(NULL, w->fl, k, 1, t2, t1, QOP_ODD, QOP_EVEN);
  QOP_wilson_diaginv_qdp(NULL, w->fl, k, t1, t2, QOP_ODD);
  QOP_wilson_dslash_qdp(NULL, w->fl, k, 1, t2, t1, QOP_EVEN, QOP_ODD);
  QDP_D_eq_D_minus_D(qd->field, qs->field, t2, QDP_even);
  QDP_destroy_D(t1);
  QDP_destroy_D(t2);
#endif
  return 0;
#undef NC
}

static int
qopqdp_wilson_precDdag(lua_State *L)
{
  BEGIN_ARGS;
  GET_WILSON(w);
  GET_QOPQDP_DFERMION(qd);
  GET_QOPQDP_DFERMION(qs);
  GET_DOUBLE(mass);
  END_ARGS;
  wilson_set(w, 2);
  double k = kappa(w,mass);
  QOP_wilson_dslash_qdp(NULL, w->fl, k, -1, qd->field, qs->field, QOP_ODD, QOP_EVEN);
  QOP_wilson_dslash_qdp(NULL, w->fl, k, -1, qd->field, qd->field, QOP_EVEN, QOP_ODD);
  QLA_Real s = -4*k*k;
  QDP_D_eq_r_times_D_plus_D(qd->field, &s, qd->field, qs->field, QDP_even);
  return 0;
}

static int
qopqdp_wilson_solve(lua_State *L)
{
#define NC QDP_get_nc(qs->field)
  QOP_evenodd_t eo=QOP_EVENODD;
  //int prec = 1;
  int max_iter = -1;
  int restart = 500;
  int max_restarts = 5;
  int precNE = 0;
  //int dag = 0;
  BEGIN_ARGS;
  GET_WILSON(w);
  GET_QOPQDP_DFERMION(qd);
  GET_QOPQDP_DFERMION(qs);
  GET_DOUBLE(mass);
  GET_DOUBLE(resid);

  if(nargs>=nextarg && lua_isstring(L, nextarg)) {
    lua_pushvalue(L, nextarg);
    const char *s = luaL_checkstring(L, -1);
    lua_pop(L, 1);
    if(strcmp(s,"precNE")==0) {
      precNE = 1;
      eo = QOP_EVEN;
    } else if(strcmp(s,"prec")==0) {
      precNE = 2;
      eo = QOP_EVEN;
    } else {
      eo = qopqdp_check_evenodd(L, nextarg);
      //printf0("unknown solver type %s\n", s);
      //qerror(1);
    }
    nextarg++;
  }
  //printf("precNE = %i\n", precNE);

  if(nargs>=nextarg && !lua_isnil(L,nextarg)) {
    if(!lua_istable(L,nextarg)) {
      qlerror0(L,1,"expecting solver paramter table\n");
    }
#define seti(s) lua_getfield(L,nextarg,#s);if(!lua_isnil(L,-1))s=luaL_checkint(L,-1);lua_pop(L,1)
    //seti(prec);
    seti(max_iter);
    seti(restart);
    seti(max_restarts);
#undef seti
    nextarg++;
  }
  if(max_iter<0) max_iter = restart*max_restarts;

  END_ARGS;

  wilson_set(w, 2);
  QOP_invert_arg_t invarg = QOP_INVERT_ARG_DEFAULT;
  invarg.max_iter = max_iter;
  invarg.restart = restart;
  invarg.max_restarts = max_restarts;
  invarg.evenodd = eo;
  QOP_resid_arg_t resarg = QOP_RESID_ARG_DEFAULT;
  resarg.rsqmin = resid*resid;
  //QDP_D_eq_zero(qd->field, QDP_all);
#if 0
  {
    //QLA_DiracFermion(*v) = QDP_expose_V(qs->field);
    //for(int i=0; i<QDP_subset_len(QDP_all); i++) {
    int i;
    QDP_loop_sites(i, QDP_all, {
	if(i==0) {
	  printf("%4i", i);
	  for(int j=0; j<QLA_Nc; j++) {
	    QLA_DiracFermion(*v) = QDP_site_ptr_readonly_V(qs->field, i);
	    printf(" %10g", QLA_real(QLA_elem_V(*v, j)));
	    printf(" %10g", QLA_imag(QLA_elem_V(*v, j)));
	  }
	  printf("\n");
	} else {
	  QLA_V_eq_zero(QDP_site_ptr_readwrite_V(qs->field,i));
	}
      });
    //QDP_reset_V(qs->field);
  }
#endif
  QOP_info_t info;
  double k = kappa(w,mass);
  if(precNE==1) {
    if(w->mg) {
      wilson_set(w, 1);
#ifdef QOP_wilsonMgSetLinks
      QOP_F_wilsonMgSetLinks(w->mg, w->ffl);
#else
      QOP_wilsonMgSetLinks(w->mg, w->ffl);
#endif
      QDP_DiracFermion *ts = QDP_create_D();
      QDP_DiracFermion *td = QDP_create_D();
      // scale by 1/4kappa^2 and solve D
      QLA_Real s = 0.25/(k*k);
      QDP_D_eq_r_times_D(ts, &s, qs->field, QDP_even);
      QDP_D_eq_zero(ts, QDP_odd);
      QDP_D_eq_zero(td, QDP_all);
#ifndef QOP_wilsonMgSolve
#define QOP_wilsonMgSolve QOP_D3_wilsonMgSolve
#endif
      QOP_wilsonMgSolve(&info, w->mg, w->fl, &invarg, &resarg, k, td, ts);
      // solve D^dag
      QDP_D_eq_gamma_times_D(ts, td, 15, QDP_even);
      QDP_D_eq_zero(ts, QDP_odd);
      QDP_D_eq_zero(td, QDP_all);
      QOP_wilsonMgSolve(&info, w->mg, w->fl, &invarg, &resarg, k, td, ts);
      QDP_D_eq_gamma_times_D(qd->field, td, 15, QDP_even);
      QDP_destroy_D(td);
      QDP_destroy_D(ts);
    } else {
      QDP_D_eq_gamma_times_D(qs->field, qs->field, 15, QDP_even);
      QOP_D_wilson_invert_ne_qdp(&info, w->fl, &invarg, &resarg, k, qd->field, qs->field);
      QDP_D_eq_gamma_times_D(qs->field, qs->field, 15, QDP_even);
      QDP_D_eq_gamma_times_D(qd->field, qd->field, 15, QDP_even);
    }
  } else {
    if(w->mg) {
      wilson_set(w, 1);
#ifdef QOP_wilsonMgSetLinks
      QOP_F_wilsonMgSetLinks(w->mg, w->ffl);
#else
      QOP_wilsonMgSetLinks(w->mg, w->ffl);
#endif
      if(invarg.evenodd!=QOP_EVENODD) {
	QDP_DiracFermion *ts = QDP_create_D();
	QDP_DiracFermion *td = QDP_create_D();
	QDP_D_eq_zero(td, QDP_all);
	if(invarg.evenodd==QOP_EVEN) {
	  QDP_D_eq_zero(ts, QDP_odd);
	  QDP_D_eq_D(ts, qs->field, QDP_even);
	} else {
	  QDP_D_eq_zero(ts, QDP_even);
	  QDP_D_eq_D(ts, qs->field, QDP_odd);
	}
	QOP_wilsonMgSolve(&info,w->mg,w->fl,&invarg,&resarg,k,td,ts);
	if(invarg.evenodd==QOP_EVEN) {
	  QDP_D_eq_D(qd->field, td, QDP_even);
	} else {
	  QDP_D_eq_D(qd->field, td, QDP_odd);
	}
	QDP_destroy_D(td);
	QDP_destroy_D(ts);
      } else {
	QOP_wilsonMgSolve(&info,w->mg,w->fl,&invarg,&resarg,k,qd->field,qs->field);
      }
    } else {
      //wilsonInvert(&info, fla, &invarg, rap, mass, nm, qqd, qs->field);
      switch(QOP_Nc) {
#ifdef HAVE_NC3
      case 3: QOP_D3_wilson_invert_qdp(&info,(QOP_3_FermionLinksWilson*)w->fl,&invarg,&resarg,k,(QDP_3_DiracFermion*)qd->field,(QDP_3_DiracFermion*)qs->field);
	break;
#endif
      default:
	QOP_D_wilson_invert_qdp(&info, w->fl, &invarg, &resarg,k,qd->field,qs->field);
      }
    }
    if(precNE==2) {
      QLA_Real s = 0.5/k;
      QDP_D_eq_r_times_D(qd->field, &s, qd->field, QDP_even);
    }
  }
#if 0
  {
    int i;
    QDP_loop_sites(i, QDP_all, {
	QLA_Real nrm2;
	QLA_r_eq_norm2_V(&nrm2, QDP_site_ptr_readonly_V(*qqd,i));
	if(nrm2>0) {
	  printf("%4i", i);
	  int x[4];
	  QDP_get_coords(x, 0, i);
	  printf(" %i %i %i %i", x[0], x[1], x[2], x[3]);
	  for(int j=0; j<QLA_Nc; j++) {
	    QLA_DiracFermion(*v) = QDP_site_ptr_readonly_V(*qqd, i);
	    printf(" %10g", QLA_real(QLA_elem_V(*v, j)));
	    printf(" %10g", QLA_imag(QLA_elem_V(*v, j)));
	  }
	  printf("\n");
	}
      });
  }
  exit(0);
#endif
  w->time = info.final_sec;
  w->flops = info.final_flop;
  w->its = resarg.final_iter;
  w->rsq = resarg.final_rsq;
  return 0;
#undef NC
}

static int
qopqdp_wilson_force(lua_State *L)
{
  int prec = 2;
  int deriv = 0;
  int all = 0;
  BEGIN_ARGS;
  GET_WILSON(w);
  GET_GAUGE(f);
  GET_AS_QOPQDP_DFERMION_ARRAY(nql, ql);
  GET_AS_QOPQDP_DFERMION_ARRAY(nqr, qr);
  GET_DOUBLE_ARRAY(nm, ms);
  GET_DOUBLE_ARRAY(ne, eps);
  if(nextarg<=nargs) { // options table
    tableLoopKeys(L, nextarg) {
      tableLoopKeysIfKeySetInt(L, "prec", prec);
      else tableLoopKeysIfKeySetInt(L, "deriv", deriv);
      else tableLoopKeysIfKeySetInt(L, "all", all);
      tableLoopKeysEnd(L);
    }
    nextarg++;
  }
  END_ARGS;
  qassert(nql==nqr);
  qassert(nql==nm);
  qassert(nql==ne);

  QOP_info_t info;
  prec = 2; // 1 not supported yet
  if(prec==1) {
#if 0
    wilson_set(h, 1);
    // factor of 4 normalizes force to that of phi^+ [s-4*Deo*Doe]^-1 phi
    QLA_F_Real qeps[ne]; for(int i=0; i<ne; i++) qeps[i] = 4*eps[i];
    QDP_F_DiracFermion *qcv[nq];
    for(int i=0; i<nq; i++) {
      qcv[i] = QDP_F_create_D();
      QDP_FD_D_eq_D(qcv[i], q[i]->field, QDP_all);
    }
    QDP_F_ColorMatrix *fforce[f->nd];
    for(int i=0; i<f->nd; i++) {
      fforce[i] = QDP_F_create_M();
      QDP_F_M_eq_zero(fforce[i], QDP_all);
    }
    QOP_F_Force *qf = QOP_F_create_F_from_qdp(fforce);
    //QOP_F_wilson_force_multi_qdp(&info, h->ffl, qf, &h->coeffs, qeps, qcv, &nq);
    QOP_F_extract_F_to_qdp(fforce, qf);
    for(int i=0; i<f->nd; i++) {
      QDP_DF_M_eq_M(f->links[i], fforce[i], QDP_all);
      QDP_F_destroy_M(fforce[i]);
    }
    for(int i=0; i<nq; i++) {
      QDP_F_destroy_D(qcv[i]);
    }
    QOP_F_destroy_F(qf);
#endif
  } else {
    wilson_set(w, 2);
    QLA_Real qks[nql], qeps[nql];
    QDP_DiracFermion *dfl[nql], *dfr[nql];
    for(int i=0; i<nql; i++) {
      qks[i] = kappa(w,ms[i]);
      qeps[i] = eps[i];
      dfl[i] = ql[i]->field;
      dfr[i] = qr[i]->field;
    }
    //for(int i=0; i<f->nd; i++) QDP_M_eq_zero(f->force[i], QDP_all);
    if(deriv) {
      rephase_F(f);
      if(w->coeffs.aniso!=1) {
	QLA_Real a = 1/w->coeffs.aniso;
	for(int i=0; i<f->nd-1; i++) {
	  QDP_M_eq_r_times_M(f->links[i], &a, f->links[i], QDP_all);
	}
      }
      if(all) {
	QOP_wilson_deriv_multi_qdp(&info, w->fl, f->links, qeps, dfl, dfr, nql);
      } else {
	QOP_wilson_deriv_prec_multi_qdp(&info, w->fl, f->links, qks, qeps, dfl, dfr, nql);
      }
      rephase_F(f);
      if(w->coeffs.aniso!=1) {
	QLA_Real a = w->coeffs.aniso;
	for(int i=0; i<f->nd-1; i++) {
	  QDP_M_eq_r_times_M(f->links[i], &a, f->links[i], QDP_all);
	}
      }
    } else {
      if(all) {
	QOP_wilson_force_multi_qdp(&info, w->fl, f->links, qeps, dfl, dfr, nql);
      } else {
	QOP_wilson_force_prec_multi_qdp(&info, w->fl, f->links, qks, qeps, dfl, dfr, nql);
      }
    }
  }

  f->time = info.final_sec;
  f->flops = info.final_flop;
  return 0;
}

static int
qopqdp_wilson_quark(lua_State *L)
{
  lattice_t *lat = NULL;
  int nc;
  BEGIN_ARGS;
  GET_WILSON(w);
  END_ARGS;
  if(w->g) lat = w->g->lat;
  else lat = qopqdp_get_default_lattice(L);
  if(w->g) nc = w->g->nc;
  else nc = lat->defaultNc;
  if(QOP_Precision=='F') {
    qopqdp_dfermionF_create(L, nc, lat);
  } else {
    qopqdp_dfermionD_create(L, nc, lat);
  }
  return 1;
}

static int
qopqdp_wilson_time(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  wilson_t *w = qopqdp_wilson_check(L, 1);
  lua_pushnumber(L, w->time);
  return 1;
}

static int
qopqdp_wilson_flops(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  wilson_t *w = qopqdp_wilson_check(L, 1);
  lua_pushnumber(L, w->flops);
  return 1;
}

static int
qopqdp_wilson_its(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  wilson_t *w = qopqdp_wilson_check(L, 1);
  lua_pushnumber(L, w->its);
  return 1;
}

static int
qopqdp_wilson_rsq(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  wilson_t *w = qopqdp_wilson_check(L, 1);
  lua_pushnumber(L, w->rsq);
  return 1;
}

static int
qopqdp_wilson_printcoeffs(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  //wilson_t *w = qopqdp_wilson_check(L, 1);
  //if(QDP_this_node==0) {
  //printf0("wilson coeffs:\n");
  //}
  return 0;
}

static struct luaL_Reg wilson_reg[] = {
  { "__gc",     qopqdp_wilson_gc },
  { "set",      qopqdp_wilson_set },
  { "mgSet",    qopqdp_wilson_mg_set },
  { "mgSetup",  qopqdp_wilson_mg_setup },
  { "quark",    qopqdp_wilson_quark },
  { "D",        qopqdp_wilson_D },
  { "Ddag",     qopqdp_wilson_Ddag },
  { "precD",    qopqdp_wilson_precD },
  { "precDdag", qopqdp_wilson_precDdag },
  { "solve",    qopqdp_wilson_solve },
  { "force",    qopqdp_wilson_force },
  { "time",     qopqdp_wilson_time },
  { "flops",    qopqdp_wilson_flops },
  { "its",      qopqdp_wilson_its },
  { "rsq",      qopqdp_wilson_rsq },
  { "printcoeffs", qopqdp_wilson_printcoeffs },
  { NULL, NULL}
};

wilson_t *
qopqdp_wilson_create(lua_State *L, int nc, lattice_t *lat)
{
  wilson_t *w = lua_newuserdata(L, sizeof(wilson_t));
  w->time = 0;
  w->flops = 0;
  w->rsq = 0;
  w->its = 0;
  w->nc = nc;
  w->lat = lat;
  w->coeffs = QOP_WILSON_COEFFS_ZERO;
  w->fl = NULL;
  w->ffl = NULL;
  w->mg = NULL;
  w->g = NULL;
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, wilson_reg);
  }
  lua_setmetatable(L, -2);
  return w;
}
