#include <string.h>
#include "qhmc_qopqdp_common.h"

static char *mtname = "qopqdp.dw";

dw_t *
qopqdp_dw_check(lua_State *L, int idx)
{
  luaL_checkudata(L, idx, mtname);
  dw_t *w = lua_touserdata(L, idx);
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
qopqdp_dw_free(lua_State *L, int idx)
{
  dw_t *w = qopqdp_dw_check(L, idx);
  if(w->fl) {
    QOP_dw_destroy_L(w->fl);
    w->fl = NULL;
  }
  if(w->ffl) {
    QOP_F_dw_destroy_L(w->ffl);
    w->ffl = NULL;
  }
}

static int
qopqdp_dw_gc(lua_State *L)
{
  qopqdp_dw_free(L, -1);
  return 0;
}

static void
qopqdp_dw_set_coeffs(QOP_dw_coeffs_t *coeffs)
{
  // none for now
}

static void
qopqdp_dw_set_opts(void)
{
  QOP_opt_t opt[3];
  opt[0].tag = "st";
  opt[0].value = 3;
  opt[1].tag = "ns";
  opt[1].value = 8;
  opt[2].tag = "nm";
  opt[2].value = 8;
  QOP_dw_invert_set_opts(opt, 3);
  //opt[0].tag = "fnmat_src_min";
  //opt[0].value = 0;
  //QOP_dw_force_set_opts(opt, 1);
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

static void
dw_set(dw_t *w, int prec)
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
  bcphase[nd-1].re = -1;
  int signmask[4] = {0,0,0,0};
  QOP_staggered_sign_t ss;
  ss.signmask = signmask;
  int r0[4] = {0,0,0,0};
  qopqdp_dw_set_coeffs(&w->coeffs);
  QOP_info_t info;
  if(prec==1) {
    QDP_F_ColorMatrix *flinks[nd];
    for(int i=0; i<nd; i++) {
      flinks[i] = QDP_F_create_M();
      QDP_FD_M_eq_M(flinks[i], g->links[i], QDP_all);
    }
    QOP_F_GaugeField *qg = QOP_F_create_G_from_qdp(flinks);
    QOP_F_rephase_G(qg, r0, &bc, &ss);
    w->ffl = QOP_F_dw_create_L_from_G(&info, &w->coeffs, qg);
    QOP_F_destroy_G(qg);
    for(int i=0; i<nd; i++) {
      QDP_F_destroy_M(flinks[i]);
    }
  } else {
    QOP_GaugeField *qg = QOP_create_G_from_qdp(g->links);
    QOP_rephase_G(qg, r0, &bc, &ss);
    w->fl = QOP_dw_create_L_from_G(&info, &w->coeffs, qg);
    QOP_destroy_G(qg);
  }
  w->time = info.final_sec;
  w->flops = info.final_flop;
#undef NC
}

static int
qopqdp_dw_set(lua_State *L)
{
  int nargs = lua_gettop(L);
  qassert(nargs>=2 && nargs<=3);
  dw_t *w = qopqdp_dw_check(L, 1);
  gauge_t *g = qopqdp_gauge_check(L, 2);
  double prec = luaL_optint(L, 3, 0);
  if(w->fl) {
    QOP_dw_destroy_L(w->fl);
    w->fl = NULL;
  }
  if(w->ffl) {
    QOP_F_dw_destroy_L(w->ffl);
    w->ffl = NULL;
  }
  qopqdp_dw_set_opts();
  w->g = g;
  if(prec==1) dw_set(w, 1);
  if(prec==2) dw_set(w, 2);
  return 0;
}

static int
qopqdp_dw_D(lua_State *L)
{
  BEGIN_ARGS;
  GET_DW(dw);
  GET_AS_QOPQDP_DFERMION_ARRAY(nqd, qqd);
  GET_AS_QOPQDP_DFERMION_ARRAY(nqs, qqs);
  GET_DOUBLE(M5);
  GET_DOUBLE(mf);
  OPT_EVENODD(eod, QOP_EVENODD);
  OPT_EVENODD(eos, QOP_EVENODD);
  END_ARGS;
  int ls = nqd;
  qassert(ls==nqs);
  QDP_DiracFermion *qd[ls], *qs[ls];
  for(int i=0; i<ls; i++) {
    qd[i] = qqd[i]->field;
    qs[i] = qqs[i]->field;
  }
  dw_set(dw, 2);
  QOP_dw_dslash_qdp(NULL, dw->fl, M5, mf, 1, qd, qs, ls, eod, eos);
  return 0;
}

static int
qopqdp_dw_Ddag(lua_State *L)
{
  BEGIN_ARGS;
  GET_DW(dw);
  GET_AS_QOPQDP_DFERMION_ARRAY(nqd, qqd);
  GET_AS_QOPQDP_DFERMION_ARRAY(nqs, qqs);
  GET_DOUBLE(M5);
  GET_DOUBLE(mf);
  OPT_EVENODD(eod, QOP_EVENODD);
  OPT_EVENODD(eos, QOP_EVENODD);
  END_ARGS;
  int ls = nqd;
  qassert(ls==nqs);
  QDP_DiracFermion *qd[ls], *qs[ls];
  for(int i=0; i<ls; i++) {
    qd[i] = qqd[i]->field;
    qs[i] = qqs[i]->field;
  }
  dw_set(dw, 2);
  QOP_dw_dslash_qdp(NULL, dw->fl, M5, mf, -1, qd, qs, ls, eod, eos);
  return 0;
}

static int
qopqdp_dw_precD(lua_State *L)
{
  BEGIN_ARGS;
  GET_DW(dw);
  GET_AS_QOPQDP_DFERMION_ARRAY(nqd, qqd);
  GET_AS_QOPQDP_DFERMION_ARRAY(nqs, qqs);
  GET_DOUBLE(M5);
  GET_DOUBLE(mf);
  END_ARGS;
  int ls = nqd;
  qassert(ls==nqs);
  QDP_DiracFermion *qd[ls], *qs[ls];
  for(int i=0; i<ls; i++) {
    qd[i] = qqd[i]->field;
    qs[i] = qqs[i]->field;
  }
  dw_set(dw, 2);
  QOP_dw_dslash_qdp(NULL, dw->fl, M5, mf, 1, qd, qs, ls, QOP_ODD, QOP_EVEN);
  QOP_dw_dslash_qdp(NULL, dw->fl, M5, mf, 1, qd, qd, ls, QOP_EVEN, QOP_ODD);
  //QDP_D_eq_r_times_D_plus_D(qd->df, &s, qd->df, qs->df, QDP_even);
  return 0;
}

static int
qopqdp_dw_precDdag(lua_State *L)
{
  BEGIN_ARGS;
  GET_DW(dw);
  GET_AS_QOPQDP_DFERMION_ARRAY(nqd, qqd);
  GET_AS_QOPQDP_DFERMION_ARRAY(nqs, qqs);
  GET_DOUBLE(M5);
  GET_DOUBLE(mf);
  END_ARGS;
  int ls = nqd;
  qassert(ls==nqs);
  QDP_DiracFermion *qd[ls], *qs[ls];
  for(int i=0; i<ls; i++) {
    qd[i] = qqd[i]->field;
    qs[i] = qqs[i]->field;
  }
  dw_set(dw, 2);
  QOP_dw_dslash_qdp(NULL, dw->fl, M5, mf, -1, qd, qs, ls, QOP_ODD, QOP_EVEN);
  QOP_dw_dslash_qdp(NULL, dw->fl, M5, mf, -1, qd, qd, ls, QOP_EVEN, QOP_ODD);
  //QDP_D_eq_r_times_D_plus_D(qd->df, &s, qd->df, qs->df, QDP_even);
  return 0;
}

static int
qopqdp_dw_solve(lua_State *L)
{
  BEGIN_ARGS;
  GET_DW(w);
  GET_AS_QOPQDP_DFERMION_ARRAY(nqd, qqd);
  GET_AS_QOPQDP_DFERMION_ARRAY(nqs, qqs);
  GET_DOUBLE(M5);
  GET_DOUBLE(mf);
  GET_DOUBLE(resid);
  OPT_STRING(sprec, NULL);
  OPT_TABLE_LEN_INDEX(nopt, iopt);
  END_ARGS;
  int ls = nqd;
  qassert(ls==nqs);
  QDP_DiracFermion *qd[ls], *qs[ls];
  for(int i=0; i<ls; i++) {
    qd[i] = qqd[i]->field;
    qs[i] = qqs[i]->field;
  }

  QOP_evenodd_t eo=QOP_EVENODD;
  //int prec = 1;
  int max_iter = -1;
  int restart = 500;
  int max_restarts = 5;
  int precNE = 0;
  //int dag = 0;
  if(sprec) {
    if(strcmp(sprec,"precNE")==0) {
      precNE = 1;
      eo = QOP_EVEN;
    } else if(strcmp(sprec,"prec")==0) {
      precNE = 2;
      eo = QOP_EVEN;
    } else {
      eo = qopqdp_check_evenodd(L, 7);
    }
  }
  //printf("precNE = %i\n", precNE);
  if(iopt) {
#define seti(s) lua_getfield(L,iopt,#s);if(!lua_isnil(L,-1))s=luaL_checkint(L,-1);lua_pop(L,1)
    //seti(prec);
    seti(max_iter);
    seti(restart);
    seti(max_restarts);
#undef seti
  }
  if(max_iter<0) max_iter = restart*max_restarts;

  dw_set(w, 2);
  QOP_invert_arg_t invarg = QOP_INVERT_ARG_DEFAULT;
  invarg.max_iter = max_iter;
  invarg.restart = restart;
  invarg.max_restarts = max_restarts;
  invarg.evenodd = eo;
  QOP_resid_arg_t resarg = QOP_RESID_ARG_DEFAULT;
  resarg.rsqmin = resid*resid;
  for(int s=0; s<ls; s++) QDP_D_eq_zero(qd[s], QDP_all);
  QOP_info_t info;
  if(precNE==1) {
    //for(int s=0; s<qd->ls; s++) QDP_D_eq_gamma_times_D(qs->df[s], qs->df[s], 15, QDP_even);
    //QOP_D_dw_invert_ne_qdp(&info, w->fl, &invarg, &resarg, mass, qd->df, qs->df);
    //for(int s=0; s<qd->ls; s++) QDP_D_eq_gamma_times_D(qs->df[s], qs->df[s], 15, QDP_even);
    //for(int s=0; s<qd->ls; s++) QDP_D_eq_gamma_times_D(qd->df[s], qd->df[s], 15, QDP_even);
  } else if(precNE==2) {
    QOP_D_dw_invert_qdp(&info, w->fl, &invarg, &resarg, M5, mf, qd, qs, ls);
    //QLA_Real f = 0.5/k;
    //for(int s=0; s<qd->ls; s++) QDP_D_eq_r_times_D(qd->df[s], &f, qd->df[s], QDP_even);
  } else {
    //dwInvert(&info, fla, &invarg, rap, mass, nm, qqd, qs->df);
    QOP_D_dw_invert_qdp(&info, w->fl, &invarg, &resarg, M5, mf, qd, qs, ls);
  }
  w->time = info.final_sec;
  w->flops = info.final_flop;
  w->its = resarg.final_iter;
  return 0;
}

static int
qopqdp_dw_force(lua_State *L)
{
#if 0 // FIXME
  BEGIN_ARGS;
  GET_DW(w);
  GET_GAUGE(f);
  GET_TABLE_LEN_INDEX(nql, iql);
  GET_TABLE_LEN_INDEX(nqr, iqr);
  qassert(nql==nqr);
  //GET_DOUBLE(M5);
  GET_AS_DOUBLE_ARRAY(nm, ms);
  GET_AS_DOUBLE_ARRAY(ne, eps);
  //OPT_INT(prec, 1);
  END_ARGS;
  int ls = nqd;
  qassert(ls==nqs);
  QDP_DiracFermion *qd[ls], *qs[ls];
  for(int i=0; i<ls; i++) {
    qd[i] = qqd[i]->field;
    qs[i] = qqs[i]->field;
  }

  int narg = lua_gettop(L);
  qassert(narg>=6 && narg<=7);
  dw_t *w = qopqdp_dw_check(L, 1);
  gauge_t *f = qopqdp_gauge_check(L, 2);
  int nql; get_table_len(L, 3, &nql);
  dwquark_t *ql[nql]; qopqdp_dwquark_array_check(L, 3, nql, ql);
  int nqr; get_table_len(L, 4, &nqr);
  dwquark_t *qr[nqr]; qopqdp_dwquark_array_check(L, 4, nqr, qr);

  int nm; get_table_len(L, 5, &nm);
  qassert(nql==nm);
  double ms[nm]; qhmc_get_double_array(L, 5, nm, ms);

  int ne; get_table_len(L, 6, &ne);
  qassert(nql==ne);
  double eps[ne]; qhmc_get_double_array(L, 6, ne, eps);
  int prec = luaL_optint(L, 7, 1);

  QOP_info_t info;
  if(prec==1) {
#if 0
    dw_set(h, 1);
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
    //QOP_F_dw_force_multi_qdp(&info, h->ffl, qf, &h->coeffs, qeps, qcv, &nq);
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
    dw_set(w, 2);
    //QLA_Real qks[nql], qeps[nql];
    //QDP_DiracFermion *dfl[nql], *dfr[nql];
    //for(int i=0; i<nql; i++) {
      //qks[i] = kappa(ms[i]);
      //qeps[i] = eps[i];
      //dfl[i] = ql[i]->df[0];
      //dfr[i] = qr[i]->df[0];
    //}
    for(int i=0; i<f->nd; i++) QDP_M_eq_zero(f->links[i], QDP_all);
    //QOP_dw_force_prec_multi_qdp(&info, w->fl, f->force, qks, qeps, dfl, dfr, nql);
  }

  f->time = info.final_sec;
  f->flops = info.final_flop;
#endif
  return 0;
}

static int
qopqdp_dw_quark(lua_State *L)
{
  lattice_t *lat = NULL;
  int nc;
  BEGIN_ARGS;
  GET_DW(dw);
  END_ARGS;
  if(dw->g) lat = dw->g->lat;
  else lat = qopqdp_get_default_lattice(L);
  if(dw->g) nc = dw->g->nc;
  else nc = lat->defaultNc;
  lua_newtable(L);
  for(int i=0; i<dw->ls; i++) {
    lua_pushnumber(L, i+1);
    if(QOP_Precision=='F') {
      qopqdp_dfermionF_create(L, nc, lat);
    } else {
      qopqdp_dfermionD_create(L, nc, lat);
    }
    lua_settable(L, -3);
  }
  return 1;
}

static int
qopqdp_dw_time(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  dw_t *w = qopqdp_dw_check(L, 1);
  lua_pushnumber(L, w->time);
  return 1;
}

static int
qopqdp_dw_flops(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  dw_t *w = qopqdp_dw_check(L, 1);
  lua_pushnumber(L, w->flops);
  return 1;
}

static int
qopqdp_dw_its(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  dw_t *w = qopqdp_dw_check(L, 1);
  lua_pushnumber(L, w->its);
  return 1;
}

static int
qopqdp_dw_printcoeffs(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  //dw_t *w = qopqdp_dw_check(L, 1);
  QOP_dw_coeffs_t coeffs;
  qopqdp_dw_set_coeffs(&coeffs);
  //if(QDP_this_node==0) {
  //printf0("dw coeffs:\n");
  //}
  return 0;
}

static struct luaL_Reg dw_reg[] = {
  { "__gc",     qopqdp_dw_gc },
  { "set",      qopqdp_dw_set },
  { "quark",    qopqdp_dw_quark },
  { "D",        qopqdp_dw_D },
  { "Ddag",     qopqdp_dw_Ddag },
  { "precD",    qopqdp_dw_precD },
  { "precDdag", qopqdp_dw_precDdag },
  { "solve",    qopqdp_dw_solve },
  { "force",    qopqdp_dw_force },
  { "time",     qopqdp_dw_time },
  { "flops",    qopqdp_dw_flops },
  { "its",      qopqdp_dw_its },
  { "printcoeffs", qopqdp_dw_printcoeffs },
  { NULL, NULL}
};

dw_t *
qopqdp_dw_create(lua_State *L, int nc, lattice_t *lat)
{
  //if(lat==NULL) lat = qopqdp_get_default_lattice(L);
  //if(nc==0) nc = lat->defaultNc;
  int ls = luaL_checkinteger(L, 1);
  dw_t *w = lua_newuserdata(L, sizeof(dw_t));
  w->ls = ls;
  w->fl = NULL;
  w->ffl = NULL;
  w->g = NULL;
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, dw_reg);
  }
  lua_setmetatable(L, -2);
  return w;
}
