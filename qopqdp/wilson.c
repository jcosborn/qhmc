#include <string.h>
#include "qhmc_qopqdp_common.h"

static char *mtname = "qopqdp.wilson";

#define kappa(m) (0.5/(4+m))

wilson_t *
qopqdp_wilson_check(lua_State *L, int idx)
{
  luaL_checkudata(L, idx, mtname);
  wilson_t *w = lua_touserdata(L, idx);
#if 0
  int hasmt = lua_getmetatable(L, idx);
  qassert(hasmt==1);
  luaL_getmetatable(L, mtname);
  int eq = lua_equal(L, -1, -2);
  qassert(eq==1);
  lua_pop(L, 2);
#endif
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
qopqdp_wilson_set_coeffs(QOP_wilson_coeffs_t *coeffs)
{
  coeffs->clov_s = 0;
  coeffs->clov_t = 0;
  coeffs->aniso = 1;
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
  //opt[0].tag = "fnmat_src_min";
  //opt[0].value = 0;
  //QOP_wilson_force_set_opts(opt, 1);
#ifdef __bg__
  QDP_set_block_size(64);
#else
  QDP_set_block_size(1024);
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

static void
wilson_set(wilson_t *w, int prec)
{
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
  bcphase[nd-1].re = -1;
  int signmask[4] = {0,0,0,0};
  QOP_staggered_sign_t ss;
  ss.signmask = signmask;
  int r0[4] = {0,0,0,0};
  qopqdp_wilson_set_coeffs(&w->coeffs);
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
}

static int
qopqdp_wilson_set(lua_State *L)
{
  int nargs = lua_gettop(L);
  qassert(nargs>=2 && nargs<=3);
  wilson_t *w = qopqdp_wilson_check(L, 1);
  gauge_t *g = qopqdp_gauge_check(L, 2);
  double prec = luaL_optint(L, 3, 0);
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

static int
qopqdp_wilson_D(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==4 || narg==6);
  wilson_t *w = qopqdp_wilson_check(L, 1);
  wquark_t *qd = qopqdp_wquark_check(L, 2);
  wquark_t *qs = qopqdp_wquark_check(L, 3);
  double mass = luaL_checknumber(L, 4);
  QOP_evenodd_t eod=QOP_EVENODD, eos=QOP_EVENODD;
  if(narg!=4) {
    eod = qopqdp_check_evenodd(L, 5);
    eos = qopqdp_check_evenodd(L, 6);
  }
  wilson_set(w, 2);
  QOP_wilson_dslash_qdp(NULL, w->fl, kappa(mass), 1, qd->df, qs->df, eod, eos);
  return 0;
}

static int
qopqdp_wilson_Ddag(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==4 || narg==6);
  wilson_t *w = qopqdp_wilson_check(L, 1);
  wquark_t *qd = qopqdp_wquark_check(L, 2);
  wquark_t *qs = qopqdp_wquark_check(L, 3);
  double mass = luaL_checknumber(L, 4);
  QOP_evenodd_t eod=QOP_EVENODD, eos=QOP_EVENODD;
  if(narg!=4) {
    eod = qopqdp_check_evenodd(L, 5);
    eos = qopqdp_check_evenodd(L, 6);
  }
  wilson_set(w, 2);
  QOP_wilson_dslash_qdp(NULL, w->fl, kappa(mass), -1, qd->df, qs->df, eod, eos);
  return 0;
}

static int
qopqdp_wilson_precD(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==4);
  wilson_t *w = qopqdp_wilson_check(L, 1);
  wquark_t *qd = qopqdp_wquark_check(L, 2);
  wquark_t *qs = qopqdp_wquark_check(L, 3);
  double mass = luaL_checknumber(L, 4);
  wilson_set(w, 2);
  double k = kappa(mass);
  QOP_wilson_dslash_qdp(NULL, w->fl, k, 1, qd->df, qs->df, QOP_ODD, QOP_EVEN);
  QOP_wilson_dslash_qdp(NULL, w->fl, k, 1, qd->df, qd->df, QOP_EVEN, QOP_ODD);
  QLA_Real s = -4*k*k;
  QDP_D_eq_r_times_D_plus_D(qd->df, &s, qd->df, qs->df, QDP_even);
  return 0;
}

static int
qopqdp_wilson_precDdag(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==4);
  wilson_t *w = qopqdp_wilson_check(L, 1);
  wquark_t *qd = qopqdp_wquark_check(L, 2);
  wquark_t *qs = qopqdp_wquark_check(L, 3);
  double mass = luaL_checknumber(L, 4);
  wilson_set(w, 2);
  double k = kappa(mass);
  QOP_wilson_dslash_qdp(NULL, w->fl, k, -1, qd->df, qs->df, QOP_ODD, QOP_EVEN);
  QOP_wilson_dslash_qdp(NULL, w->fl, k, -1, qd->df, qd->df, QOP_EVEN, QOP_ODD);
  QLA_Real s = -4*k*k;
  QDP_D_eq_r_times_D_plus_D(qd->df, &s, qd->df, qs->df, QDP_even);
  return 0;
}

static int
qopqdp_wilson_solve(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg>=5 && narg<=7);
  wilson_t *w = qopqdp_wilson_check(L, 1);
  wquark_t *qd = qopqdp_wquark_check(L, 2);
  wquark_t *qs = qopqdp_wquark_check(L, 3);
  double mass = luaL_checknumber(L, 4);
  double resid = luaL_checknumber(L, 5);
  QOP_evenodd_t eo=QOP_EVENODD;
  //int prec = 1;
  int max_iter = -1;
  int restart = 500;
  int max_restarts = 5;
  int nextarg = 6;
  int precNE = 0;
  //int dag = 0;
  if(narg>=nextarg && lua_isstring(L, nextarg)) {
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
  if(narg>=nextarg && !lua_isnil(L,nextarg)) {
    if(!lua_istable(L,nextarg)) {
      printf0("expecting solver paramter table\n");
      qerror(1);
    }
#define seti(s) lua_getfield(L,nextarg,#s);if(!lua_isnil(L,-1))s=luaL_checkint(L,-1);lua_pop(L,1)
    //seti(prec);
    seti(max_iter);
    seti(restart);
    seti(max_restarts);
#undef seti
  }
  if(max_iter<0) max_iter = restart*max_restarts;

  wilson_set(w, 2);
  QOP_invert_arg_t invarg = QOP_INVERT_ARG_DEFAULT;
  invarg.max_iter = max_iter;
  invarg.restart = restart;
  invarg.max_restarts = max_restarts;
  invarg.evenodd = eo;
  QOP_resid_arg_t resarg = QOP_RESID_ARG_DEFAULT;
  resarg.rsqmin = resid*resid;
  QDP_D_eq_zero(qd->df, QDP_all);
#if 0
  {
    //QLA_DiracFermion *v = QDP_expose_V(qs->df);
    //for(int i=0; i<QDP_subset_len(QDP_all); i++) {
    int i;
    QDP_loop_sites(i, QDP_all, {
	if(i==0) {
	  printf("%4i", i);
	  for(int j=0; j<QLA_Nc; j++) {
	    QLA_DiracFermion *v = QDP_site_ptr_readonly_V(qs->df, i);
	    printf(" %10g", QLA_real(QLA_elem_V(*v, j)));
	    printf(" %10g", QLA_imag(QLA_elem_V(*v, j)));
	  }
	  printf("\n");
	} else {
	  QLA_V_eq_zero(QDP_site_ptr_readwrite_V(qs->df,i));
	}
      });
    //QDP_reset_V(qs->df);
  }
#endif
  QOP_info_t info;
  double k = kappa(mass);
  if(precNE==1) {
    QDP_D_eq_gamma_times_D(qs->df, qs->df, 15, QDP_even);
    QOP_D_wilson_invert_ne_qdp(&info, w->fl, &invarg, &resarg, k, qd->df, qs->df);
    QDP_D_eq_gamma_times_D(qs->df, qs->df, 15, QDP_even);
    QDP_D_eq_gamma_times_D(qd->df, qd->df, 15, QDP_even);
  } else if(precNE==2) {
    QOP_D_wilson_invert_qdp(&info, w->fl, &invarg, &resarg, k, qd->df, qs->df);
    QLA_Real s = 0.5/k;
    QDP_D_eq_r_times_D(qd->df, &s, qd->df, QDP_even);
  } else {
    //wilsonInvert(&info, fla, &invarg, rap, mass, nm, qqd, qs->df);
    QOP_D_wilson_invert_qdp(&info, w->fl, &invarg, &resarg, k, qd->df, qs->df);
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
	    QLA_DiracFermion *v = QDP_site_ptr_readonly_V(*qqd, i);
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
  return 0;
}

static int
qopqdp_wilson_force(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg>=6 && narg<=7);
  wilson_t *w = qopqdp_wilson_check(L, 1);
  force_t *f = qopqdp_force_check(L, 2);
  int nql; get_table_len(L, 3, &nql);
  wquark_t *ql[nql]; qopqdp_wquark_array_check(L, 3, nql, ql);
  int nqr; get_table_len(L, 4, &nqr);
  qassert(nql==nqr);
  wquark_t *qr[nqr]; qopqdp_wquark_array_check(L, 4, nqr, qr);
  int nm; get_table_len(L, 5, &nm);
  qassert(nql==nm);
  double ms[nm]; get_double_array(L, 5, nm, ms);
  int ne; get_table_len(L, 6, &ne);
  qassert(nql==ne);
  double eps[ne]; get_double_array(L, 6, ne, eps);
  int prec = luaL_optint(L, 7, 1);

  QOP_info_t info;
  if(prec==1) {
#if 0
    wilson_set(h, 1);
    // factor of 4 normalizes force to that of phi^+ [s-4*Deo*Doe]^-1 phi
    QLA_F_Real qeps[ne]; for(int i=0; i<ne; i++) qeps[i] = 4*eps[i];
    QDP_F_DiracFermion *qcv[nq];
    for(int i=0; i<nq; i++) {
      qcv[i] = QDP_F_create_D();
      QDP_FD_D_eq_D(qcv[i], q[i]->df, QDP_all);
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
      QDP_DF_M_eq_M(f->force[i], fforce[i], QDP_all);
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
      qks[i] = kappa(ms[i]);
      qeps[i] = eps[i];
      dfl[i] = ql[i]->df;
      dfr[i] = qr[i]->df;
    }
    for(int i=0; i<f->nd; i++) QDP_M_eq_zero(f->force[i], QDP_all);
    QOP_wilson_force_prec_multi_qdp(&info, w->fl, f->force, qks, qeps, dfl, dfr, nql);
  }

  f->time = info.final_sec;
  f->flops = info.final_flop;
  return 0;
}

static int
qopqdp_wilson_quark(lua_State* L)
{
  qassert(lua_gettop(L)==1);
  qopqdp_wquark_create(L);
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
qopqdp_wilson_printcoeffs(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  //wilson_t *w = qopqdp_wilson_check(L, 1);
  QOP_wilson_coeffs_t coeffs;
  qopqdp_wilson_set_coeffs(&coeffs);
  //if(QDP_this_node==0) {
  //printf0("wilson coeffs:\n");
  //}
  return 0;
}

static struct luaL_Reg wilson_reg[] = {
  { "__gc",     qopqdp_wilson_gc },
  { "set",      qopqdp_wilson_set },
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
  { "printcoeffs", qopqdp_wilson_printcoeffs },
  { NULL, NULL}
};

wilson_t *
qopqdp_wilson_create(lua_State* L)
{
  wilson_t *w = lua_newuserdata(L, sizeof(wilson_t));
  w->fl = NULL;
  w->ffl = NULL;
  w->g = NULL;
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, wilson_reg);
  }
  lua_setmetatable(L, -2);
  return w;
}
