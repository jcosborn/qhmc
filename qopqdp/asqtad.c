#include <string.h>
#include "qhmc_qopqdp_common.h"

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
  opt[0].value = 1;
  opt[1].tag = "ns";
  opt[1].value = 8;
  opt[2].tag = "nm";
  opt[2].value = 16;
  QOP_asqtad_invert_set_opts(opt, 3);
  opt[0].tag = "fnmat_src_min";
  opt[0].value = 0;
  QOP_asqtad_force_set_opts(opt, 1);
#ifdef __bg__
  QDP_set_block_size(64);
#else
  QDP_set_block_size(1024);
#endif
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
  if((prec==1&&h->ffl)||(prec==2&&h->fl)) {
    return;
  }
  gauge_t *g = h->g;
  int nd = g->nd;
  QOP_info_t info;
  if(prec==1) {
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
  } else {
    QOP_GaugeField *qg = QOP_create_G_from_qdp(g->links);
    QOP_rephase_G(qg, h->r0, &h->bc, &h->ssign);
    h->fl = QOP_asqtad_create_L_from_G(&info, &h->coeffs, qg);
    QOP_destroy_G(qg);
  }
  h->time = info.final_sec;
  h->flops = info.final_flop;
}

static int
qopqdp_asqtad_set(lua_State *L)
{
  int nargs = lua_gettop(L);
  qassert(nargs>=2&&nargs<=3);
  asqtad_t *h = qopqdp_asqtad_check(L, 1);
  gauge_t *g = qopqdp_gauge_check(L, 2);
  double prec = luaL_optint(L, 3, 0);
  if(h->fl) {
    QOP_asqtad_destroy_L(h->fl);
    h->fl = NULL;
  }
  if(h->ffl) {
    QOP_F_asqtad_destroy_L(h->ffl);
    h->ffl = NULL;
  }
  qopqdp_asqtad_set_opts();
  h->g = g;
  if(prec==1) asqtad_set(h, 1);
  if(prec==2) asqtad_set(h, 2);
  return 0;
}

static int
qopqdp_asqtad_D(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==4 || narg==6);
  asqtad_t *h = qopqdp_asqtad_check(L, 1);
  squark_t *qd = qopqdp_squark_check(L, 2);
  squark_t *qs = qopqdp_squark_check(L, 3);
  double mass = luaL_checknumber(L, 4);
  QOP_evenodd_t eod=QOP_EVENODD, eos=QOP_EVENODD;
  if(narg!=4) {
    eod = qopqdp_check_evenodd(L, 5);
    eos = qopqdp_check_evenodd(L, 6);
  }
  asqtad_set(h, 2);
  QOP_asqtad_dslash_qdp(NULL, h->fl, mass, qd->cv, qs->cv, eod, eos);
  return 0;
}

static int
qopqdp_asqtad_Ddag(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==4 || narg==6);
  asqtad_t *h = qopqdp_asqtad_check(L, 1);
  squark_t *qd = qopqdp_squark_check(L, 2);
  squark_t *qs = qopqdp_squark_check(L, 3);
  double mass = luaL_checknumber(L, 4);
  QOP_evenodd_t eod=QOP_EVENODD, eos=QOP_EVENODD;
  if(narg!=4) {
    eod = qopqdp_check_evenodd(L, 5);
    eos = qopqdp_check_evenodd(L, 6);
  }
  asqtad_set(h, 2);
  QOP_asqtad_dslash_qdp(NULL, h->fl, -mass, qd->cv, qs->cv, eod, eos);
  QDP_V_eqm_V(qd->cv, qd->cv, QDP_all);
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
  int narg = lua_gettop(L);
  qassert(narg>=5 && narg<=7);
  asqtad_t *h = qopqdp_asqtad_check(L, 1);
  int nqd; get_table_len(L, 2, &nqd);
  squark_t *qd[nqd]; qopqdp_squark_array_check(L, 2, nqd, qd);
  squark_t *qs = qopqdp_squark_check(L, 3);
  int nm; get_table_len(L, 4, &nm);
  qassert(nqd==nm);
  double mass[nm]; get_double_array(L, 4, nm, mass);
  double resid = luaL_checknumber(L, 5);
  QOP_evenodd_t eo=QOP_EVENODD;
  //int prec = 1;
  int max_iter = -1;
  int restart = 500;
  int max_restarts = 5;
  int use_prev_soln = 0;
  int nextarg = 6;
  if(narg>=nextarg && lua_isstring(L, nextarg)) {
    eo = qopqdp_check_evenodd(L, nextarg);
    nextarg++;
  }
  if(narg>=nextarg && !lua_isnil(L,nextarg)) {
    if(!lua_istable(L,nextarg)) {
      qerror0("expecting solver paramter table\n");
    }
#define seti(s) lua_getfield(L,nextarg,#s);if(!lua_isnil(L,-1))s=luaL_checkint(L,-1);lua_pop(L,1)
    //seti(prec);
    seti(max_iter);
    seti(restart);
    seti(max_restarts);
    seti(use_prev_soln);
#undef seti
  }
  if(max_iter<0) max_iter = restart*max_restarts;

  asqtad_set(h, 2);
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
    qqd[i] = qd[i]->cv;
    if(!use_prev_soln) {
      QDP_V_eq_zero(qqd[i], QDP_all);
    }
  }
#if 0
  {
    //QLA_ColorVector *v = QDP_expose_V(qs->cv);
    //for(int i=0; i<QDP_subset_len(QDP_all); i++) {
    int i;
    QDP_loop_sites(i, QDP_all, {
	if(i==0) {
	  printf("%4i", i);
	  for(int j=0; j<QLA_Nc; j++) {
	    QLA_ColorVector *v = QDP_site_ptr_readonly_V(qs->cv, i);
	    printf(" %10g", QLA_real(QLA_elem_V(*v, j)));
	    printf(" %10g", QLA_imag(QLA_elem_V(*v, j)));
	  }
	  printf("\n");
	} else {
	  QLA_V_eq_zero(QDP_site_ptr_readwrite_V(qs->cv,i));
	}
      });
    //QDP_reset_V(qs->cv);
  }
#endif
  QOP_info_t info;
  asqtadInvert(&info, h->fl, &invarg, rap, mass, nm, qqd, qs->cv);
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
	    QLA_ColorVector *v = QDP_site_ptr_readonly_V(*qqd, i);
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
  return 0;
}

static int
qopqdp_asqtad_force(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg>=4 && narg<=5);
  asqtad_t *h = qopqdp_asqtad_check(L, 1);
  force_t *f = qopqdp_force_check(L, 2);
  int nq; get_table_len(L, 3, &nq);
  squark_t *q[nq]; qopqdp_squark_array_check(L, 3, nq, q);
  int ne; get_table_len(L, 4, &ne);
  qassert(nq==ne);
  double eps[ne]; get_double_array(L, 4, ne, eps);
  int prec = 1;
  int deriv = 0;
  if(narg==5) { // options table
    tableLoopKeys(L, 5) {
      tableLoopKeysIfKeySetInt(L, "prec", prec);
      else tableLoopKeysIfKeySetInt(L, "deriv", deriv);
      tableLoopKeysEnd(L);
    }
  }

  QOP_info_t info;
  if(prec==1) {
    asqtad_set(h, 1);
    // factor of 4 normalizes force to that of phi^+ [s-4*Deo*Doe]^-1 phi
    QLA_F_Real qeps[ne]; for(int i=0; i<ne; i++) qeps[i] = 4*eps[i];
    QDP_F_ColorVector *qcv[nq];
    for(int i=0; i<nq; i++) {
      qcv[i] = QDP_F_create_V();
      QDP_FD_V_eq_V(qcv[i], q[i]->cv, QDP_all);
    }
    QDP_F_ColorMatrix *fforce[f->nd];
    for(int i=0; i<f->nd; i++) {
      fforce[i] = QDP_F_create_M();
      QDP_F_M_eq_zero(fforce[i], QDP_all);
    }
    QDP_F_ColorMatrix *flinks[f->nd];
    for(int i=0; i<f->nd; i++) {
      flinks[i] = QDP_F_create_M();
      QDP_FD_M_eq_M(flinks[i], h->g->links[i], QDP_all);
    }
    QOP_F_rephase_G_qdp(flinks, h->r0, &h->bc, &h->ssign);
    if(deriv) {
      QOP_F_asqtad_deriv_multi_qdp(&info, flinks, fforce, &h->coeffs, qeps, qcv, nq);
      QOP_F_rephase_G_qdp(fforce, h->r0, &h->bc, &h->ssign);
    } else {
      QOP_F_asqtad_force_multi_qdp(&info, flinks, fforce, &h->coeffs, qeps, qcv, nq);
    }
    for(int i=0; i<f->nd; i++) {
      QDP_DF_M_eq_M(f->force[i], fforce[i], QDP_all);
      QDP_F_destroy_M(fforce[i]);
      QDP_F_destroy_M(flinks[i]);
    }
    for(int i=0; i<nq; i++) {
      QDP_F_destroy_V(qcv[i]);
    }
  } else {
    asqtad_set(h, 2);
    // factor of 4 normalizes force to that of phi^+ [s-4*Deo*Doe]^-1 phi
    QLA_Real qeps[ne]; for(int i=0; i<ne; i++) qeps[i] = 4*eps[i];
    QDP_ColorVector *qcv[nq]; for(int i=0; i<nq; i++) qcv[i] = q[i]->cv;
    for(int i=0; i<f->nd; i++) QDP_M_eq_zero(f->force[i], QDP_all);
    QOP_rephase_G_qdp(h->g->links, h->r0, &h->bc, &h->ssign);
    if(deriv) {
      QOP_asqtad_deriv_multi_qdp(&info, h->g->links, f->force, &h->coeffs, qeps, qcv, nq);
      QOP_rephase_G_qdp(f->force, h->r0, &h->bc, &h->ssign);
    } else {
      QOP_asqtad_force_multi_qdp(&info, h->g->links, f->force, &h->coeffs, qeps, qcv, nq);
    }
    QOP_rephase_G_qdp(h->g->links, h->r0, &h->bc, &h->ssign);
  }

  f->time = info.final_sec;
  f->flops = info.final_flop;
  return 0;
}

static int
qopqdp_asqtad_squark(lua_State* L)
{
  qassert(lua_gettop(L)==1);
  qopqdp_squark_create(L);
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

static struct luaL_Reg asqtad_reg[] = {
  { "__gc",        qopqdp_asqtad_gc },
  { "set",         qopqdp_asqtad_set },
  { "quark",       qopqdp_asqtad_squark },
  { "D",           qopqdp_asqtad_D },
  { "Ddag",        qopqdp_asqtad_Ddag },
  { "solve",       qopqdp_asqtad_solve },
  { "force",       qopqdp_asqtad_force },
  { "time",        qopqdp_asqtad_time },
  { "flops",       qopqdp_asqtad_flops },
  { "its",         qopqdp_asqtad_its },
  { "coeffs",      qopqdp_asqtad_coeffs },
  { NULL, NULL}
};

asqtad_t *
qopqdp_asqtad_create(lua_State* L)
{
  asqtad_t *h = lua_newuserdata(L, sizeof(asqtad_t));
  h->fl = NULL;
  h->ffl = NULL;
  h->g = NULL;
  h->fatphase = NULL;
  h->longphase = NULL;
  h->nd = QDP_ndim();
  qopqdp_origin_alloc(&h->r0, h->nd);
  qopqdp_origin_default(h->r0, h->nd);
  qopqdp_bcphase_alloc(&h->bc, h->nd);
  qopqdp_bcphase_default(&h->bc, h->nd);
  qopqdp_signmask_alloc(&h->ssign, h->nd);
  qopqdp_signmask_default(&h->ssign, h->nd);
  qopqdp_asqtad_coeffs_default(&h->coeffs, 1);
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, asqtad_reg);
  }
  lua_setmetatable(L, -2);
  return h;
}
