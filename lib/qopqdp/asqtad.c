//#include <omp.h>
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

static char *mtname = "qopqdp.asqtad";

// bcphase

static void
qopqdp_origin_alloc(int **r0, int nd)
{
  *r0 = malloc(nd*sizeof(int));
}

static void
qopqdp_origin_free(int *r0)
{
  free(r0);
}

static void
qopqdp_origin_default(int *r0, int nd)
{
  for(int i=0; i<nd; i++) r0[i] = 0;
}

static void
qopqdp_bcphase_alloc(QOP_bc_t *bc, int nd)
{
  bc->phase = malloc(nd*sizeof(QOP_Complex));
}

static void
qopqdp_bcphase_free(QOP_bc_t *bc)
{
  if(bc->phase) free(bc->phase);
}

static void
qopqdp_bcphase_default(QOP_bc_t *bc, int nd)
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
static int
qopqdp_bcphase_check(lua_State *L, int idx, QOP_bc_t *c)
{
  if(lua_type(L,idx)!=LUA_TTABLE) {
    return 0;
  }
  return 1;
}

static int
qopqdp_bcphase_push(lua_State *L, QOP_bc_t *c, char *prefix, int maketable)
{
  if(maketable) {
    lua_createtable(L, 6, 0);
  }
  return 1;
}
#endif

static int
qopqdp_asqtad_bc(lua_State *L)
{
  BEGIN_ARGS;
  GET_ASQTAD(a);
  OPT_AS_COMPLEX_ARRAY(nbc,bc,0,NULL);
  END_ARGS;
  int nd = a->lat->nd;
  if(nbc>0) { // set
    int n = nd<nbc ? nd : nbc;
    for(int i=0; i<n; i++) {
      a->bc.phase[i].re = bc[i].r;
      a->bc.phase[i].im = bc[i].i;
    }
    return 0;
  }
  // else get
  qhmc_complex_t b[nd];
  for(int i=0; i<nd; i++) {
    b[i].r = a->bc.phase[i].re;
    b[i].i = a->bc.phase[i].im;
  }
  qhmc_push_complex_array(L, nd, b);
  return 1;
}

static void
qopqdp_signmask_alloc(QOP_staggered_sign_t *s, int nd)
{
  s->signmask = malloc(nd*sizeof(int));
}

static void
qopqdp_signmask_free(QOP_staggered_sign_t *s)
{
  free(s->signmask);
}

static void
qopqdp_signmask_default(QOP_staggered_sign_t *s, int nd)
{
#ifdef QHMC_REPRO_UNIFORM
  for(int i=0; i<nd; i++) s->signmask[i] = 0;
#else
  //int signmask[4] = {8,9,11,0};
  int k = 1;
  for(int i=0; i<nd-1; i++) {
    int t = k - (k&15) + ((k&1)<<3) + ((k&14)>>1);  // move 0 bit to 3rd
    s->signmask[i] = t;
    //printf("i: %i  k: %i  t: %i\n", i, k, t);
    k = 2*k + 1;
  }
  s->signmask[nd-1] = 0;
  //for(int i=0; i<nd; i++) {
  //printf("signmask[%i] = %i\n", i, s->signmask[i]);
  //}
#endif
}

// asqtad

asqtad_t *
qopqdp_asqtad_check(lua_State *L, int idx)
{
  luaL_checkudata(L, idx, mtname);
  asqtad_t *h = lua_touserdata(L, idx);
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
qopqdp_asqtad_free(lua_State *L, int idx)
{
  asqtad_t *h = qopqdp_asqtad_check(L, idx);
  if(h->fl) {
    QOP_asqtad_destroy_L(h->fl);
    h->fl = NULL;
  }
  if(h->ffl) {
    QOP_F_asqtad_destroy_L(h->ffl);
    h->ffl = NULL;
  }
  qopqdp_origin_free(h->r0);
  qopqdp_bcphase_free(&h->bc);
  qopqdp_signmask_free(&h->ssign);
}

static int
qopqdp_asqtad_gc(lua_State *L)
{
  qopqdp_asqtad_free(L, -1);
  return 0;
}

// asqtad coeffs

static void
qopqdp_asqtad_coeffs_zero(QOP_asqtad_coeffs_t *coeffs)
{
  coeffs->one_link = 0;
  coeffs->three_staple = 0;
  coeffs->five_staple = 0;
  coeffs->seven_staple = 0;
  coeffs->lepage = 0;
  coeffs->naik = 0;
}

static void
qopqdp_asqtad_coeffs_default(QOP_asqtad_coeffs_t *coeffs, double u0)
{
  double u2, u4;
  u2 = 1.0/(u0*u0);
  u4 = u2*u2;
#ifdef QHMC_REPRO_UNIFORM
  u2 = -u2;
#endif
  coeffs->one_link = 5.0/8.0;
  coeffs->three_staple = -u2/16.0;
  coeffs->five_staple = u4/64.0;
  coeffs->seven_staple = -(u4*u2)/384.0;
  coeffs->lepage = -u4/16.0;
  coeffs->naik = -u2/24.0;
}

static int
qopqdp_asqtad_coeffs_check(lua_State *L, int idx, QOP_asqtad_coeffs_t *c)
{
  if(lua_type(L,idx)!=LUA_TTABLE) {
    return 0;
  }
  tableLoopKeys(L, idx)	{
#define set(x) tableLoopKeysIfKeySetDouble(L, #x, c->x)
    set(one_link);
    else set(three_staple);
    else set(five_staple);
    else set(seven_staple);
    else set(lepage);
    else set(naik);
    tableLoopKeysEnd(L);
#undef set
  }
  return 1;
}

static int
qopqdp_asqtad_coeffs_push(lua_State *L, QOP_asqtad_coeffs_t *c, char *prefix, int maketable)
{
  if(maketable) {
    lua_createtable(L, 6, 0);
  }
#define push(x) lua_pushnumber(L, c->x); lua_setfield(L, -2, #x)
  push(one_link);
  push(three_staple);
  push(five_staple);
  push(seven_staple);
  push(lepage);
  push(naik);
#undef push
  return 1;
}

// asqtad set

static void
qopqdp_asqtad_set_opts(void)
{
  QOP_opt_t opt[3];
  opt[0].tag = "st";
  opt[1].tag = "ns";
  opt[2].tag = "nm";
#ifdef __bg__
  opt[0].value = 0;
  opt[1].value = 8;
  opt[2].value = 8;
  //QDP_set_block_size(1024);
#else
  opt[0].value = 1;
  opt[1].value = 8;
  opt[2].value = 16;
  //QDP_set_block_size(1024);
#endif
  QOP_asqtad_invert_set_opts(opt, 3);
  opt[0].tag = "fnmat_src_min";
  opt[0].value = 0;
  QOP_asqtad_force_set_opts(opt, 1);
  {
    QOP_opt_t opt = {.tag="reunit_allow_svd",.value=0};
    QOP_hisq_links_set_opts(&opt, 1);
  }
  {
    QOP_opt_t opt = {.tag="force_filter",.value=0};
    QOP_hisq_force_set_opts(&opt, 1);
  }
}

static void
asqtad_set(asqtad_t *h, int prec)
{
#define NC QDP_get_nc(g->links[0])
  if((prec==1&&h->ffl)||(prec==2&&h->fl)) {
    return;
  }
  gauge_t *g = h->g;
  gauge_t *g2 = h->gLong;
  QOP_info_t info;
  if(prec==1) {
#if 0
    int nd = g->nd;
    QDP_F_ColorMatrix *flinks[nd];
    for(int i=0; i<nd; i++) {
      flinks[i] = QDP_F_create_M();
      QDP_FD_M_eq_M(flinks[i], g->links[i], QDP_all);
    }
    QOP_F_GaugeField *qg = QOP_F_create_G_from_qdp(flinks);
    QOP_F_rephase_G(qg, h->r0, &h->bc, &h->ssign);
    h->ffl = QOP_F_asqtad_create_L_from_G(&info, &h->coeffs, qg);
    QOP_F_destroy_G(qg);
    for(int i=0; i<nd; i++) {
      QDP_F_destroy_M(flinks[i]);
    }
#else
    QOP_F_GaugeField *qg = QOP_FD_create_G_from_qdp(g->links);
    QOP_F_rephase_G(qg, h->r0, &h->bc, &h->ssign);
    h->ffl = QOP_F_asqtad_create_L_from_G(&info, &h->coeffs, qg);
    QOP_F_destroy_G(qg);
#endif
  } else {
    QOP_GaugeField *qg = QOP_create_G_from_qdp(g->links);
    QOP_rephase_G(qg, h->r0, &h->bc, &h->ssign);
    QOP_GaugeField *qg2 = qg;
    if(g2!=g) {
      qg2 = QOP_create_G_from_qdp(g2->links);
      QOP_rephase_G(qg2, h->r0, &h->bc, &h->ssign);
    }
    h->fl = QOP_asqtad_create_L_from_G2(&info, &h->coeffs, qg, qg2);
    QOP_destroy_G(qg);
    if(g2!=g) QOP_destroy_G(qg2);
  }
  h->time = info.final_sec;
  h->flops = info.final_flop;
#undef NC
}

static int
qopqdp_asqtad_set(lua_State *L)
{
  BEGIN_ARGS;
  GET_ASQTAD(a);
  GET_GAUGE(g);
  OPT_GAUGE(g2, g);
  OPT_INT(prec, 0);
  END_ARGS;
  if(a->fl) {
    QOP_asqtad_destroy_L(a->fl);
    a->fl = NULL;
  }
  if(a->ffl) {
    QOP_F_asqtad_destroy_L(a->ffl);
    a->ffl = NULL;
  }
  qopqdp_asqtad_set_opts();
  a->g = g;
  a->gLong = g2;
  if(prec==1) asqtad_set(a, 1);
  if(prec==2) asqtad_set(a, 2);
  return 0;
}

static int
qopqdp_asqtad_D(lua_State *L)
{
  BEGIN_ARGS;
  GET_ASQTAD(a);
  GET_QOPQDP_CVECTOR(qd);
  GET_QOPQDP_CVECTOR(qs);
  GET_DOUBLE(mass);
  OPT_EVENODD(eod, QOP_EVENODD);
  OPT_EVENODD(eos, QOP_EVENODD);
  END_ARGS;
  asqtad_set(a, 2);
  //printf("%i %i %i\n", h->g->nc, QDP_get_nc(qd->field), QDP_get_nc(qs->field));
  QOP_asqtad_dslash_qdp(NULL, a->fl, mass, qd->field, qs->field, eod, eos);
  return 0;
}

static int
qopqdp_asqtad_Ddag(lua_State *L)
{
  BEGIN_ARGS;
  GET_ASQTAD(a);
  GET_QOPQDP_CVECTOR(qd);
  GET_QOPQDP_CVECTOR(qs);
  GET_DOUBLE(mass);
  OPT_EVENODD(eod, QOP_EVENODD);
  OPT_EVENODD(eos, QOP_EVENODD);
  END_ARGS;
  asqtad_set(a, 2);
  QOP_asqtad_dslash_qdp(NULL, a->fl, -mass, qd->field, qs->field, eod, eos);
  QDP_V_eqm_V(qd->field, qd->field, QDP_all_L(qd->qlat));
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
qopqdp_asqtad_solve(lua_State *L)
{
#define NC QDP_get_nc(qs->field)
  BEGIN_ARGS;
  GET_ASQTAD(a);
  GET_AS_QOPQDP_CVECTOR_ARRAY(nqd, qd);
  GET_QOPQDP_CVECTOR(qs);
  GET_DOUBLE_ARRAY(nm, mass);
  qassert(nqd==nm);
  GET_DOUBLE(resid);

  //int ompnt = omp_get_max_threads();
  //omp_set_num_threads(1);

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

  asqtad_set(a, 2);
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
      QDP_V_eq_zero(qqd[i], QDP_all_L(qd[i]->qlat));
    }
  }
#if 0
  {
    //QLA_ColorVector(*v) = QDP_expose_V(qs->field);
    //for(int i=0; i<QDP_subset_len(QDP_all); i++) {
    int i;
    QDP_loop_sites(i, QDP_all_L(qs->qlat), {
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
    case 1: QOP_1_asqtad_solve_multi_qdp(&info, (QOP_1_FermionLinksAsqtad*)a->fl,&invarg,rap,mass,(QDP_1_ColorVector**)qqd,(QDP_1_ColorVector**)srcs,nm);
      break;
#endif
#ifdef HAVE_NC2
    case 2: QOP_2_asqtad_solve_multi_qdp(&info, (QOP_2_FermionLinksAsqtad*)a->fl,&invarg,rap,mass,(QDP_2_ColorVector**)qqd,(QDP_2_ColorVector**)srcs,nm);
      break;
#endif
#ifdef HAVE_NC3
      case 3: QOP_3_asqtad_solve_multi_qdp(&info, (QOP_3_FermionLinksAsqtad*)a->fl,&invarg,rap,mass,(QDP_3_ColorVector**)qqd,(QDP_3_ColorVector**)srcs,nm);
      break;
#endif
    default:
      QOP_asqtad_solve_multi_qdp(&info,a->fl,&invarg,rap,mass,qqd,srcs,nm);
    }
  } else
#endif
    {
      //printf0("calling asqtadInvert\n");
      asqtadInvert(&info, a->fl, &invarg, rap, mass, nm, qqd, qs->field);
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
  if(a->fatphase || a->longphase) {
    int nd = QDP_ndim();
    QDP_Complex *fp[nd], *lp[nd];
    for(int i=0; i<nd; i++) {
      fp[i] = NULL;
      if(a->fatphase[i]) {
	fp[i] = QDP_create_C();
	C_eq_1_div_C(fp[i], a->fatphase[i]);
      }
      lp[i] = NULL;
      if(a->longphase[i]) {
	lp[i] = QDP_create_C();
	C_eq_1_div_C(lp[i], a->longphase[i]);
      }
    }
    //QOP_asqtad_rephase_field_L_qdp(fla, fp, lp);
  }
  a->time = info.final_sec;
  a->flops = info.final_flop;
  a->its = resarg[0].final_iter;
  a->rsq = resarg[0].final_rsq;
  //omp_set_num_threads(ompnt);
  return 0;
#undef NC
}

static int
qopqdp_asqtad_force(lua_State *L)
{
#define NC QDP_get_nc(f->links[0])
  BEGIN_ARGS;
  GET_ASQTAD(a);
  GET_GAUGE(f);
  GET_AS_QOPQDP_CVECTOR_ARRAY(nq, q);
  GET_DOUBLE_ARRAY(ne, eps);
  qassert(nq==ne);
  int prec = 1;
  int deriv = 0;
  if(nargs==5) { // options table
    tableLoopKeys(L, 5) {
      tableLoopKeysIfKeySetInt(L, "prec", prec);
      else tableLoopKeysIfKeySetInt(L, "deriv", deriv);
      tableLoopKeysEnd(L);
    }
    nextarg++;
  }
  END_ARGS;

  QOP_info_t info;
  if(prec==1) {
    asqtad_set(a, 1);
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
    QDP_F_ColorMatrix *flinks[f->nd];
    for(int i=0; i<f->nd; i++) {
      flinks[i] = QDP_F_create_M();
      QDP_FD_M_eq_M(flinks[i], a->g->links[i], QDP_all);
    }
    QOP_F_rephase_G_qdp(flinks, a->r0, &a->bc, &a->ssign);
    if(deriv) {
      switch(QOP_Nc) {
#ifdef HAVE_NC3
      case 3:
	QOP_F3_asqtad_deriv_multi_qdp(&info, (QDP_F3_ColorMatrix**)flinks, (QDP_F3_ColorMatrix**)fforce, &a->coeffs, qeps, (QDP_F3_ColorVector**)qcv, nq);
	break;
#endif
      default: QOP_F_asqtad_deriv_multi_qdp(&info, flinks, fforce, &a->coeffs, qeps, qcv, nq);
      }
      QOP_F_rephase_G_qdp(fforce, a->r0, &a->bc, &a->ssign);
    } else {
      QOP_F_asqtad_force_multi_qdp(&info, flinks, fforce, &a->coeffs, qeps, qcv, nq);
    }
    for(int i=0; i<f->nd; i++) {
      QDP_DF_M_eq_M(f->links[i], fforce[i], QDP_all);
      QDP_F_destroy_M(fforce[i]);
      QDP_F_destroy_M(flinks[i]);
    }
    for(int i=0; i<nq; i++) {
      QDP_F_destroy_V(qcv[i]);
    }
  } else {
    asqtad_set(a, 2);
    // factor of 4 normalizes force to that of phi^+ [s-4*Deo*Doe]^-1 phi
    QLA_Real qeps[ne]; for(int i=0; i<ne; i++) qeps[i] = 4*eps[i];
    QDP_ColorVector *qcv[nq]; for(int i=0; i<nq; i++) qcv[i] = q[i]->field;
    for(int i=0; i<f->nd; i++) QDP_M_eq_zero(f->links[i], QDP_all);
    QOP_rephase_G_qdp(a->g->links, a->r0, &a->bc, &a->ssign);
    if(deriv) {
      switch(QOP_Nc) {
#ifdef HAVE_NC3
      case 3:
	QOP_3_asqtad_deriv_multi_qdp(&info, (QDP_3_ColorMatrix**)a->g->links, (QDP_3_ColorMatrix**)f->links, &a->coeffs, qeps, (QDP_3_ColorVector**)qcv, nq);
	break;
#endif
      default: QOP_asqtad_deriv_multi_qdp(&info, a->g->links, f->links, &a->coeffs, qeps, qcv, nq);
      }
      QOP_rephase_G_qdp(f->links, a->r0, &a->bc, &a->ssign);
    } else {
      QOP_asqtad_force_multi_qdp(&info, a->g->links, f->links, &a->coeffs, qeps, qcv, nq);
    }
    QOP_rephase_G_qdp(a->g->links, a->r0, &a->bc, &a->ssign);
  }

  f->time = info.final_sec;
  f->flops = info.final_flop;
  return 0;
#undef NC
}

static int
qopqdp_asqtad_squark(lua_State *L)
{
  BEGIN_ARGS;
  GET_ASQTAD(a);
  END_ARGS;
  qopqdp_cvectorD_create(L, a->nc, a->lat);
  return 1;
}

static int
qopqdp_asqtad_time(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  asqtad_t *h = qopqdp_asqtad_check(L, 1);
  lua_pushnumber(L, h->time);
  return 1;
}

static int
qopqdp_asqtad_flops(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  asqtad_t *h = qopqdp_asqtad_check(L, 1);
  lua_pushnumber(L, h->flops);
  return 1;
}

static int
qopqdp_asqtad_its(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  asqtad_t *h = qopqdp_asqtad_check(L, 1);
  lua_pushnumber(L, h->its);
  return 1;
}

static int
qopqdp_asqtad_rsq(lua_State *L)
{
  BEGIN_ARGS;
  GET_ASQTAD(a);
  END_ARGS;
  lua_pushnumber(L, a->rsq);
  return 1;
}

static int
qopqdp_asqtad_coeffs(lua_State *L)
{
  int nargs = lua_gettop(L);
  qassert(nargs>=1&&nargs<=2);
  asqtad_t *h = qopqdp_asqtad_check(L, 1);
  if(nargs==1) { // return coefficients
    qopqdp_asqtad_coeffs_push(L, &h->coeffs, "", 1);
    return 1;
  } else { // nargs==2; set coefficients
    qopqdp_asqtad_coeffs_zero(&h->coeffs);
    int set = qopqdp_asqtad_coeffs_check(L, 2, &h->coeffs);
    qassert(set!=0);
  }
  return 0;
}

#ifdef SOLVE2
#include "solve2.c"
#endif

static struct luaL_Reg asqtad_reg[] = {
  { "__gc",        qopqdp_asqtad_gc },
  { "bc",          qopqdp_asqtad_bc },
  { "set",         qopqdp_asqtad_set },
  { "quark",       qopqdp_asqtad_squark },
  { "D",           qopqdp_asqtad_D },
  { "Ddag",        qopqdp_asqtad_Ddag },
  { "solve",       qopqdp_asqtad_solve },
#ifdef SOLVE2
  { "solve2",      qopqdp_asqtad_solve2 },
#endif
  { "force",       qopqdp_asqtad_force },
  { "time",        qopqdp_asqtad_time },
  { "flops",       qopqdp_asqtad_flops },
  { "its",         qopqdp_asqtad_its },
  { "rsq",         qopqdp_asqtad_rsq },
  { "coeffs",      qopqdp_asqtad_coeffs },
  { NULL, NULL}
};

asqtad_t *
qopqdp_asqtad_create(lua_State *L, int nc, lattice_t *lat)
{
  if(lat==NULL) lat = qopqdp_get_default_lattice(L);
  if(nc==0) nc = lat->defaultNc;
  asqtad_t *h = lua_newuserdata(L, sizeof(asqtad_t));
  h->nc = nc;
  h->lat = lat;
  h->fl = NULL;
  h->ffl = NULL;
  h->g = NULL;
  h->gLong = NULL;
  h->fatphase = NULL;
  h->longphase = NULL;
  int nd = lat->nd;
  qopqdp_origin_alloc(&h->r0, nd);
  qopqdp_origin_default(h->r0, nd);
  qopqdp_bcphase_alloc(&h->bc, nd);
  qopqdp_bcphase_default(&h->bc, nd);
  qopqdp_signmask_alloc(&h->ssign, nd);
  qopqdp_signmask_default(&h->ssign, nd);
  qopqdp_asqtad_coeffs_default(&h->coeffs, 1);
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, asqtad_reg);
  }
  lua_setmetatable(L, -2);
  return h;
}
