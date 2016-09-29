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

static char *mtname = "qopqdp.hisq";

hisq_t *
qopqdp_hisq_check(lua_State *L, int idx)
{
  luaL_checkudata(L, idx, mtname);
  hisq_t *h = lua_touserdata(L, idx);
#if 0
  int hasmt = lua_getmetatable(L, idx);
  qassert(hasmt==1);
  luaL_getmetatable(L, mtname);
  int eq = lua_equal(L, -1, -2);
  qassert(eq==1);
  lua_pop(L, 2);
#endif
  return h;
}

static void
qopqdp_hisq_free(lua_State *L, int idx)
{
  hisq_t *h = qopqdp_hisq_check(L, idx);
  if(h->fl) {
    QOP_hisq_destroy_L(h->fl);
    h->fl = NULL;
  }
  if(h->ffl) {
    QOP_F_hisq_destroy_L(h->ffl);
    h->ffl = NULL;
  }
}

static int
qopqdp_hisq_gc(lua_State *L)
{
  qopqdp_hisq_free(L, -1);
  return 0;
}

static void
qopqdp_hisq_set_coeffs(QOP_hisq_coeffs_t *coeffs, double u0, double f7lf)
{
  coeffs->n_naiks = 1;
  coeffs->eps_naik[0] = 0;
  coeffs->ugroup = QOP_UNITARIZE_U3;
  //coeffs->ugroup = QOP_UNITARIZE_SU3;
  coeffs->umethod = QOP_UNITARIZE_RATIONAL;
  //coeffs->umethod = QOP_UNITARIZE_ANALYTIC;

  double u2, u4;
  u2 = 1.0/(u0*u0);
  u4 = u2*u2;
#ifdef QHMC_REPRO_UNIFORM
  u2 = -u2;
#endif
  coeffs->fat7_one_link = (1.0+3.0*f7lf)/8.0;
  coeffs->fat7_three_staple = -u2/16.0;
  coeffs->fat7_five_staple = u4/64.0;
  coeffs->fat7_seven_staple = -(u2*u4)/384.0;
  coeffs->fat7_lepage = -u4*(f7lf)/16.0;
  coeffs->asqtad_one_link = (8.0-3.0*f7lf)/8.0;
  coeffs->asqtad_three_staple = -1.0/16.0;
  coeffs->asqtad_five_staple = 1.0/64.0;
  coeffs->asqtad_seven_staple = -1.0/384.0;
  coeffs->asqtad_lepage = -(2.0-f7lf)/16.0;
  coeffs->asqtad_naik = -1.0/24.0;
  coeffs->difference_one_link = 0;
  coeffs->difference_naik = 0;
}

static void
qopqdp_hisq_set_opts(void)
{
  QOP_opt_t opt[3];
  opt[0].tag = "st";
  opt[0].value = 1;
  opt[1].tag = "ns";
  opt[1].value = 8;
  opt[2].tag = "nm";
  opt[2].value = 16;
  QOP_asqtad_invert_set_opts(opt, 3);
  //opt[0].tag = "cg";
  //opt[0].value = 0;
  //QOP_asqtad_invert_set_opts(opt, 1);
  opt[0].tag = "fnmat_src_min";
  opt[0].value = 0;
  QOP_hisq_force_set_opts(opt, 1);
#ifdef __bg__
  //QDP_set_block_size(64);
#else
  //QDP_set_block_size(1024);
#endif
  opt[0].tag = "use_fat7_lepage";
  opt[0].value = 1;
  QOP_hisq_links_set_opts(opt, 1);
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
  //int signmask[4] = {0,1,3,7};
  int signmask[4] = {8,9,11,0};
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
hisq_set(hisq_t *h, int prec)
{
#define NC QDP_get_nc(g->links[0])
  if((prec==1&&h->ffl)||(prec==2&&h->fl)) {
    return;
  }
  gauge_t *g = h->g;
  int nd = g->nd;
  QOP_Complex bcphase[nd];
  QOP_bc_t bc;
  bc.phase = bcphase;
  for(int i=0; i<nd; i++) {
    bcphase[i].re = 1;
    bcphase[i].im = 0;
  }
#ifdef QHMC_REPRO_UNIFORM
  int signmask[4] = {0,0,0,0};
#else
  bcphase[nd-1].re = -1;
  int signmask[4] = {8,9,11,0};
#endif
  QOP_staggered_sign_t ss;
  ss.signmask = signmask;
  int r0[4] = {0,0,0,0};
  qopqdp_hisq_set_coeffs(&h->coeffs, h->u0, h->f7lf);
  QOP_info_t info;
  if(prec==1) {
    QDP_F_ColorMatrix *flinks[nd];
    for(int i=0; i<nd; i++) {
      flinks[i] = QDP_F_create_M();
      QDP_FD_M_eq_M(flinks[i], g->links[i], QDP_all);
    }
    QOP_F_GaugeField *qg = QOP_F_create_G_from_qdp(flinks);
    QOP_F_rephase_G(qg, r0, &bc, &ss);
    h->ffl = QOP_F_hisq_create_L_from_G(&info, &h->coeffs, qg);
    QOP_F_destroy_G(qg);
    for(int i=0; i<nd; i++) {
      QDP_F_destroy_M(flinks[i]);
    }
  } else {
    QOP_GaugeField *qg = QOP_create_G_from_qdp(g->links);
    QOP_rephase_G(qg, r0, &bc, &ss);
    h->fl = QOP_hisq_create_L_from_G(&info, &h->coeffs, qg);
    QOP_destroy_G(qg);
  }
  h->time = info.final_sec;
  h->flops = info.final_flop;
#undef NC
}

static int
qopqdp_hisq_set(lua_State *L)
{
  int nargs = lua_gettop(L);
  qassert(nargs>=4&&nargs<=5);
  hisq_t *h = qopqdp_hisq_check(L, 1);
  gauge_t *g = qopqdp_gauge_check(L, 2);
  double u0 = luaL_checknumber(L, 3);
  double f7lf = luaL_checknumber(L, 4);
  double prec = luaL_optint(L, 5, 0);
  if(h->fl) {
    QOP_hisq_destroy_L(h->fl);
    h->fl = NULL;
  }
  if(h->ffl) {
    QOP_F_hisq_destroy_L(h->ffl);
    h->ffl = NULL;
  }
  qopqdp_hisq_set_opts();
  h->g = g;
  h->u0 = u0;
  h->f7lf = f7lf;
  if(prec==1) hisq_set(h, 1);
  if(prec==2) hisq_set(h, 2);
  return 0;
}

static int
qopqdp_hisq_D(lua_State *L)
{
  BEGIN_ARGS;
  GET_HISQ(h);
  GET_QOPQDP_CVECTOR(qd);
  GET_QOPQDP_CVECTOR(qs);
  GET_DOUBLE(mass);
  OPT_EVENODD(eod, QOP_EVENODD);
  OPT_EVENODD(eos, QOP_EVENODD);
  END_ARGS;
  hisq_set(h, 2);
  QOP_FermionLinksAsqtad *fla = QOP_get_asqtad_links_from_hisq(h->fl)[0];
  QOP_asqtad_dslash_qdp(NULL, fla, mass, qd->field, qs->field, eod, eos);
  return 0;
}

static int
qopqdp_hisq_Ddag(lua_State *L)
{
  BEGIN_ARGS;
  GET_HISQ(h);
  GET_QOPQDP_CVECTOR(qd);
  GET_QOPQDP_CVECTOR(qs);
  GET_DOUBLE(mass);
  OPT_EVENODD(eod, QOP_EVENODD);
  OPT_EVENODD(eos, QOP_EVENODD);
  END_ARGS;
  hisq_set(h, 2);
  QOP_FermionLinksAsqtad *fla = QOP_get_asqtad_links_from_hisq(h->fl)[0];
  QOP_asqtad_dslash_qdp(NULL, fla, -mass, qd->field, qs->field, eod, eos);
  QDP_V_eqm_V(qd->field, qd->field, QDP_all);
  return 0;
}

static void
C_eq_1_div_C(QDP_Complex *r, QDP_Complex *a)
{
  int i;
  QDP_loop_sites(i, QDP_all, {
      QLA_c_eq_r_div_c(*(QLA_Complex*)QDP_site_ptr_readwrite_C(r,i), 1,
		       *(QLA_Complex*)QDP_site_ptr_readonly_C(a,i));
    });
}

static int
qopqdp_hisq_solve(lua_State *L)
{
#define NC QDP_get_nc(qs->field)
  BEGIN_ARGS;
  GET_HISQ(h);
  GET_AS_QOPQDP_CVECTOR_ARRAY(nqd, qd);
  GET_QOPQDP_CVECTOR(qs);
  GET_DOUBLE_ARRAY(nm, mass);
  qassert(nqd==nm);
  GET_DOUBLE(resid);

  QOP_evenodd_t eo=QOP_EVENODD;
  //int prec = 1;
  int max_iter = -1;
  int restart = 500;
  int max_restarts = 5;
  int use_prev_soln = 0;
  double mixed_rsq = 0;

  if(nargs>=nextarg && lua_isstring(L, nextarg)) {
    eo = qopqdp_check_evenodd(L, nextarg);
    nextarg++;
  }

  if(nargs>=nextarg && !lua_isnil(L,nextarg)) {
    if(!lua_istable(L,nextarg)) {
      qlerror0(L,1,"expecting solver paramter table\n");
    }
#define seti(s) lua_getfield(L,nextarg,#s);if(!lua_isnil(L,-1))s=luaL_checkint(L,-1);lua_pop(L,1)
#define setd(s) lua_getfield(L,nextarg,#s);if(!lua_isnil(L,-1))s=luaL_checknumber(L,-1);lua_pop(L,1)
    //seti(prec);
    seti(max_iter);
    seti(restart);
    seti(max_restarts);
    seti(use_prev_soln);
    setd(mixed_rsq);
#undef seti
#undef setd
    nextarg++;
  }
  if(max_iter<0) max_iter = restart*max_restarts;

  END_ARGS;

  hisq_set(h, 2);
  QOP_FermionLinksAsqtad *fla = QOP_get_asqtad_links_from_hisq(h->fl)[0];
  //if(h->fatphase || h->longphase) {
  //QOP_asqtad_rephase_field_L_qdp(fla, h->fatphase, h->longphase);
  //}
  QOP_invert_arg_t invarg = QOP_INVERT_ARG_DEFAULT;
  invarg.max_iter = max_iter;
  invarg.restart = restart;
  invarg.max_restarts = max_restarts;
  invarg.evenodd = eo;
  QOP_resid_arg_t resarg[nm], *rap[nm];
  for(int i=0; i<nm; i++) {
    resarg[i] = QOP_RESID_ARG_DEFAULT;
    resarg[i].rsqmin = resid*resid;
    rap[i] = &resarg[i];
  }
  QDP_ColorVector *qqd[nqd];
  for(int i=0; i<nqd; i++) {
    qqd[i] = qd[i]->field;
    if(!use_prev_soln) {
      QDP_V_eq_zero(qqd[i], QDP_all);
    }
  }
#if 0
  {
    //QLA_ColorVector(*v) = QDP_expose_V(qs->field);
    //for(int i=0; i<QDP_subset_len(QDP_all); i++) {
    int i;
    QDP_loop_sites(i, QDP_all, {
	if(i==0) {
	  printf("%4i", i);
	  for(int j=0; j<QLA_Nc; j++) {
	    QLA_ColorVector(*v) = QDP_site_ptr_readonly_V(qs->field, i);
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
#ifdef QOP_asqtad_solve_multi_qdp
  if(eo!=QOP_EVENODD) {
    invarg.mixed_rsq = mixed_rsq;
    QDP_ColorVector *srcs[nm];
    for(int i=0; i<nm; i++) srcs[i] = qs->field;
    //QOP_verbose(QOP_VERB_MED);
    switch(QOP_Nc) {
#ifdef HAVE_NC1
    case 1: QOP_1_asqtad_solve_multi_qdp(&info, (QOP_1_FermionLinksAsqtad*)fla,&invarg,rap,mass,(QDP_1_ColorVector**)qqd,(QDP_1_ColorVector**)srcs,nm);
      break;
#endif
#ifdef HAVE_NC2
    case 2: QOP_2_asqtad_solve_multi_qdp(&info, (QOP_2_FermionLinksAsqtad*)fla,&invarg,rap,mass,(QDP_2_ColorVector**)qqd,(QDP_2_ColorVector**)srcs,nm);
      break;
#endif
#ifdef HAVE_NC3
      case 3: QOP_3_asqtad_solve_multi_qdp(&info, (QOP_3_FermionLinksAsqtad*)fla,&invarg,rap,mass,(QDP_3_ColorVector**)qqd,(QDP_3_ColorVector**)srcs,nm);
      break;
#endif
    default:
      QOP_asqtad_solve_multi_qdp(&info,fla,&invarg,rap,mass,qqd,srcs,nm);
    }
  } else
#endif
    {
      asqtadInvert(&info, fla, &invarg, rap, mass, nm, qqd, qs->field);
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
	    QLA_ColorVector(*v) = QDP_site_ptr_readonly_V(*qqd, i);
	    printf(" %10g", QLA_real(QLA_elem_V(*v, j)));
	    printf(" %10g", QLA_imag(QLA_elem_V(*v, j)));
	  }
	  printf("\n");
	}
      });
  }
  exit(0);
#endif
  if(h->fatphase || h->longphase) {
    int nd = QDP_ndim();
    QDP_Complex *fp[nd], *lp[nd];
    for(int i=0; i<nd; i++) {
      fp[i] = NULL;
      if(h->fatphase[i]) {
	fp[i] = QDP_create_C();
	C_eq_1_div_C(fp[i], h->fatphase[i]);
      }
      lp[i] = NULL;
      if(h->longphase[i]) {
	lp[i] = QDP_create_C();
	C_eq_1_div_C(lp[i], h->longphase[i]);
      }
    }
    //QOP_asqtad_rephase_field_L_qdp(fla, fp, lp);
  }
  h->time = info.final_sec;
  h->flops = info.final_flop;
  h->its = resarg[0].final_iter;
  h->rsq = resarg[0].final_rsq;
  return 0;
#undef NC
}

static int
qopqdp_hisq_force(lua_State *L)
{
#define NC QDP_get_nc(f->links[0])
  BEGIN_ARGS;
  GET_HISQ(h);
  GET_GAUGE(f);
  GET_AS_QOPQDP_CVECTOR_ARRAY(nq, q);
  GET_DOUBLE_ARRAY(ne, eps);
  OPT_INT(prec, 1);
  END_ARGS;

  QOP_info_t info;
  if(prec==1) {
    hisq_set(h, 1);
    // factor of 4 normalizes force to that of phi^+ [s-4*Deo*Doe]^-1 phi
    QLA_F_Real qeps[ne]; for(int i=0; i<ne; i++) qeps[i] = 4*eps[i];
    QDP_F_ColorVector *qcv[nq];
    for(int i=0; i<nq; i++) {
      qcv[i] = QDP_F_create_V();
      QDP_FD_V_eq_V(qcv[i], q[i]->field, QDP_all);
    }
    QDP_F_ColorMatrix *fforce[f->nd];
    for(int i=0; i<f->nd; i++) {
      fforce[i] = QDP_F_create_M();
      QDP_F_M_eq_zero(fforce[i], QDP_all);
    }
    switch(QOP_Nc) {
#ifdef HAVE_NC3
    case 3: QOP_F3_hisq_force_multi_qdp(&info, (QOP_F3_FermionLinksHisq*)h->ffl, (QDP_F3_ColorMatrix**)fforce, &h->coeffs, qeps, (QDP_F3_ColorVector**)qcv, &nq);
      break;
#endif
    default: QOP_F_hisq_force_multi_qdp(&info, h->ffl, fforce, &h->coeffs, qeps, qcv, &nq);
    }
    for(int i=0; i<f->nd; i++) {
      QDP_DF_M_eq_M(f->links[i], fforce[i], QDP_all);
      QDP_F_destroy_M(fforce[i]);
    }
    for(int i=0; i<nq; i++) {
      QDP_F_destroy_V(qcv[i]);
    }
  } else {
    hisq_set(h, 2);
    // factor of 4 normalizes force to that of phi^+ [s-4*Deo*Doe]^-1 phi
    QLA_Real qeps[ne]; for(int i=0; i<ne; i++) qeps[i] = 4*eps[i];
    QDP_ColorVector *qcv[nq]; for(int i=0; i<nq; i++) qcv[i] = q[i]->field;
    for(int i=0; i<f->nd; i++) QDP_M_eq_zero(f->links[i], QDP_all);
    switch(QOP_Nc) {
#ifdef HAVE_NC3
    case 3: QOP_3_hisq_force_multi_qdp(&info, (QOP_3_FermionLinksHisq*)h->fl, (QDP_3_ColorMatrix**)f->links, &h->coeffs, qeps, (QDP_3_ColorVector**)qcv, &nq);
      break;
#endif
    default: QOP_hisq_force_multi_qdp(&info, h->fl, f->links, &h->coeffs, qeps, qcv, &nq);
    }
  }

  f->time = info.final_sec;
  f->flops = info.final_flop;
  return 0;
#undef NC
}

static int
qopqdp_hisq_squark(lua_State *L)
{
  BEGIN_ARGS;
  GET_HISQ(h);
  END_ARGS;
  qopqdp_cvector_create(L, h->nc, h->lat);
  return 1;
}

static int
qopqdp_hisq_time(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  hisq_t *h = qopqdp_hisq_check(L, 1);
  lua_pushnumber(L, h->time);
  return 1;
}

static int
qopqdp_hisq_flops(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  hisq_t *h = qopqdp_hisq_check(L, 1);
  lua_pushnumber(L, h->flops);
  return 1;
}

static int
qopqdp_hisq_its(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  hisq_t *h = qopqdp_hisq_check(L, 1);
  lua_pushnumber(L, h->its);
  return 1;
}

static int
qopqdp_hisq_rsq(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  hisq_t *h = qopqdp_hisq_check(L, 1);
  lua_pushnumber(L, h->rsq);
  return 1;
}

static int
qopqdp_hisq_printcoeffs(lua_State *L)
{
  qassert(lua_gettop(L)==3);
  //hisq_t *h = qopqdp_hisq_check(L, 1);
  double u0 = luaL_checknumber(L, 2);
  double f7lf = luaL_checknumber(L, 3);
  QOP_hisq_coeffs_t coeffs;
  qopqdp_hisq_set_coeffs(&coeffs, u0, f7lf);
  if(QDP_this_node==0) {
    printf0("hisq coeffs:\n");
#define p(s) printf0(" %-20s = % g\n", #s, coeffs.s)
    p(fat7_one_link);
    p(fat7_three_staple);
    p(fat7_five_staple);
    p(fat7_seven_staple);
    p(fat7_lepage);
    p(asqtad_one_link);
    p(asqtad_three_staple);
    p(asqtad_five_staple);
    p(asqtad_seven_staple);
    p(asqtad_lepage);
    p(asqtad_naik);
#undef p
  }
  return 0;
}

static int
qopqdp_hisq_force_filter(lua_State *L)
{
  qassert(lua_gettop(L)==2);
  //hisq_t *h = qopqdp_hisq_check(L, 1);
  double newval = luaL_checknumber(L, 2);
  QOP_opt_t opt = {.tag="force_filter",.value=newval};
  QOP_hisq_force_set_opts(&opt, 1);
  return 0;
}

static int
qopqdp_hisq_set_phases(lua_State *L)
{
  qassert(lua_gettop(L)==2);
  //hisq_t *h = qopqdp_hisq_check(L, 1);
  double newval = luaL_checknumber(L, 2);
  QOP_opt_t opt = {.tag="force_filter",.value=newval};
  QOP_hisq_force_set_opts(&opt, 1);
  return 0;
}

static struct luaL_Reg hisq_reg[] = {
  { "__gc",    qopqdp_hisq_gc },
  { "set",     qopqdp_hisq_set },
  { "quark",   qopqdp_hisq_squark },
  { "D",       qopqdp_hisq_D },
  { "Ddag",    qopqdp_hisq_Ddag },
  { "solve",   qopqdp_hisq_solve },
  { "force",   qopqdp_hisq_force },
  { "time",    qopqdp_hisq_time },
  { "flops",   qopqdp_hisq_flops },
  { "its",     qopqdp_hisq_its },
  { "rsq",     qopqdp_hisq_rsq },
  { "printcoeffs", qopqdp_hisq_printcoeffs },
  { "force_filter", qopqdp_hisq_force_filter },
  { "set_phases", qopqdp_hisq_set_phases },
  { NULL, NULL}
};

hisq_t *
qopqdp_hisq_create(lua_State *L, int nc, lattice_t *lat)
{
  if(lat==NULL) lat = qopqdp_get_default_lattice(L);
  if(nc==0) nc = lat->defaultNc;
  hisq_t *h = lua_newuserdata(L, sizeof(hisq_t));
  h->time = 0;
  h->flops = 0;
  h->rsq = 0;
  h->its = 0;
  h->nc = nc;
  h->lat = lat;
  h->fl = NULL;
  h->ffl = NULL;
  h->g = NULL;
  h->fatphase = NULL;
  h->longphase = NULL;
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, hisq_reg);
  }
  lua_setmetatable(L, -2);
  return h;
}
