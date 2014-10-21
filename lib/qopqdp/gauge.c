#include "qhmc_qopqdp_common.h"
#include <string.h>
#include <math.h>
#include <qla_d.h>

static char *gmtname = "qopqdp.gauge" STR(QOP_PrecisionLetter);

gauge_t *
qopqdp_gauge_opt(lua_State *L, int *idx, int required, gauge_t *def)
{
  gauge_t *g;
  if(required) {
    g = luaL_checkudata(L, *idx, gmtname);
    (*idx)++;
  } else {
    g = luaL_testudata(L, *idx, gmtname);
    if(g==NULL) g = def;
    else (*idx)++;
  }
  return g;
}

gauge_t *
qopqdp_gauge_check(lua_State *L, int idx)
{
  luaL_checkudata(L, idx, gmtname);
  gauge_t *g = lua_touserdata(L, idx);
  return g;
}

void
qopqdp_gauge_array_check(lua_State *L, int idx, int n, gauge_t *g[n])
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qassert(lua_objlen(L, idx)==n);
  lua_pushvalue(L, idx); // make copy for indexing convenience
  for(int i=0; i<n; i++) {
    lua_pushnumber(L, i+1);
    lua_gettable(L, -2);
    g[i] = qopqdp_gauge_check(L, -1);
    lua_pop(L, 1);
  }
  lua_pop(L, 1);
}

static void
qopqdp_gauge_free(lua_State *L, int idx)
{
  gauge_t *g = qopqdp_gauge_check(L, idx);
  for(int i=0; i<g->nd; i++) {
    QDP_destroy_M(g->links[i]);
  }
}

static int
qopqdp_gauge_gc(lua_State *L)
{
  qopqdp_gauge_free(L, -1);
  return 0;
}

static int
qopqdp_gauge_len(lua_State *L)
{
  BEGIN_ARGS;
  GET_GAUGE(g);
  nextarg++; // ignore second argument
  END_ARGS;
  lua_pushinteger(L, g->nd);
  return 1;
}

static int
qopqdp_gauge_call(lua_State *L)
{
  BEGIN_ARGS;
  GET_GAUGE(g);
  GET_INT(dim);
  OPT_QOPQDP_CMATRIX(m, NULL);
  END_ARGS;
  if(m==NULL) { // return colormatrix for direction
    qopqdp_cmatrix_wrap(L, g->lat, g->links[dim-1], 0);
    return 1;
  }
  // set direction
  QDP_M_eq_M(g->links[dim-1], m->field, QDP_all_L(g->qlat));
  return 0;
}

static int
qopqdp_gauge_nc(lua_State *L)
{
  BEGIN_ARGS;
  GET_GAUGE(g);
  END_ARGS;
  int nc = g->nc;
  lua_pushinteger(L, nc);
  return 1;
}

static int
qopqdp_gauge_lattice(lua_State *L)
{
  BEGIN_ARGS;
  GET_GAUGE(g);
  END_ARGS;
  lattice_t *lat = g->lat;
  qopqdp_lattice_wrap(L, lat->qlat, lat->defaultPrecision, lat->defaultNc, 0);
  return 1;
}

static int
qopqdp_gauge_unit(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  gauge_t *g = qopqdp_gauge_check(L, -1);
  QLA_Complex z;
  QLA_c_eq_r(z, 1);
  for(int i=0; i<g->nd; i++) {
    QDP_M_eq_c(g->links[i], &z, QDP_all_L(g->qlat));
  }
  return 0;
}

static void
gauge_random(NCPROT QLA_ColorMatrix(*m), int i, void *args)
{
  QLA_D_Complex d1, d2;
  QLA_ColorMatrix(m1);
  QLA_ColorMatrix(m2);
  QLA_ColorMatrix(m3);
  QLA_Complex c;
  QLA_RandomState *s = &((QLA_RandomState *)args)[i];
  QLA_M_eq_gaussian_S(&m1, s);
  QLA_M_eq_Ma_times_M(&m2, &m1, &m1);
  QLA_M_eq_invsqrt_M(&m3, &m2);
  if(QLA_Nc>1) {
    QLA_M_eq_M_times_M(&m2, &m1, &m3);
    QLA_C_eq_det_M(&c, &m2);
    QLA_DF_c_eq_c(d1, c);
    d2 = QLA_D_clog(&d1);
    QLA_c_eq_r_times_c(d1, -1./QLA_Nc, d2);
    d2 = QLA_D_cexp(&d1);
    QLA_FD_c_eq_c(c, d2);
    QLA_M_eq_C_times_M(m, &c, &m2);
  } else {
    QLA_M_eq_M_times_M(m, &m1, &m3);
  }
}

#if 0
static void
make_herm(NCPROT QLA_ColorMatrix(*m), int idx, void *args)
{
  QLA_Complex tr;
  QLA_c_eq_r(tr, 0);
  for(int i=0; i<QLA_Nc; i++) {
    for(int j=i; j<QLA_Nc; j++) {
      QLA_Complex t1, t2;
      QLA_c_eq_c(t1, QLA_elem_M(*m,i,j));
      QLA_c_eq_c(t2, QLA_elem_M(*m,j,i));
      QLA_c_peq_ca(t1, t2);
      QLA_c_eq_r_times_c(t1, 0.5, t1);
      QLA_c_eq_c(QLA_elem_M(*m,i,j), t1);
      QLA_c_eq_ca(QLA_elem_M(*m,j,i), t1);
    }
    QLA_c_peq_c(tr, QLA_elem_M(*m,i,i));
  }
  QLA_c_eq_r_times_c(tr, 1./QLA_Nc, tr);
  for(int i=0; i<QLA_Nc; i++) {
    QLA_c_meq_c(QLA_elem_M(*m,i,i), tr);
  }
}
#endif

static int
qopqdp_gauge_random(lua_State *L)
{
  BEGIN_ARGS;
  GET_GAUGE(g);
  END_ARGS;
  QLA_RandomState *qrs = QDP_expose_S(g->lat->rs);
  for(int i=0; i<g->nd; i++) {
    QDP_M_eq_funcia(g->links[i], gauge_random, qrs, QDP_all_L(g->qlat));
  }
  QDP_reset_S(g->lat->rs);
  return 0;
}

// 1: gauge
// 2: filename
// return file, record metadata
static int
qopqdp_gauge_load(lua_State *L)
{
#define NC nc
  BEGIN_ARGS;
  GET_GAUGE(g);
  GET_STRING(fn);
  END_ARGS;
  //printf0("loading lattice file %s\n", fn);
  //double dt = -QDP_time();
  QDP_String *fmd = QDP_string_create();
  QDP_String *rmd = QDP_string_create();
  QDP_Reader *qr = QDP_open_read_L(g->qlat, fmd, (char*)fn);
  int nd = g->nd;
  int prec, type, nc;
  qopqdp_get_prec_type_nc(qr, &prec, &type, &nc);
  if(prec=='F') {
#if QDP_Precision == 'F'
    QDP_F_vread_M(qr, rmd, g->links, nd);
#else
    QDP_F_ColorMatrix *tm[nd];
    for(int i=0; i<nd; i++) {
      tm[i] = QDP_F_create_M_L(g->qlat);
    }
    QDP_F_vread_M(qr, rmd, tm, nd);
    for(int i=0; i<nd; i++) {
      QDP_DF_M_eq_M(g->links[i], tm[i], QDP_all_L(g->qlat));
      QDP_F_destroy_M(tm[i]);
    }
#endif
  } else { // record precision is 'D'
#if QDP_Precision == 'F'
    QDP_D_ColorMatrix *tm[nd];
    for(int i=0; i<nd; i++) {
      tm[i] = QDP_D_create_M_L(g->qlat);
    }
    QDP_D_vread_M(qr, rmd, tm, nd);
    for(int i=0; i<nd; i++) {
      QDP_FD_M_eq_M(g->links[i], tm[i], QDP_all_L(g->qlat));
      QDP_D_destroy_M(tm[i]);
    }
#else
    QDP_D_vread_M(qr, rmd, g->links, nd);
#endif
  }
  QDP_close_read(qr);
  lua_pushstring(L, QDP_string_ptr(fmd));
  lua_pushstring(L, QDP_string_ptr(rmd));
  QDP_string_destroy(fmd);
  QDP_string_destroy(rmd);
  //dt += QDP_time();
  //printf0(" loaded in %g seconds\n", dt);
  return 2;
#undef NC
}

// 1: gauge
// 2: filename
// 3: file metadata
// 4: record metadata
// 5: (optional) precision
static int
qopqdp_gauge_save(lua_State *L)
{
#define NC QDP_get_nc(g->links[0])
  BEGIN_ARGS;
  GET_GAUGE(g);
  GET_STRING(fn);
  GET_STRING(fmds);
  GET_STRING(rmds);
  OPT_STRING(precision, "F");
  END_ARGS;
  //printf0("saving lattice file %s\n", fn);
  //double dt = -QDP_time();
  QDP_String *fmd = QDP_string_create();
  QDP_string_set(fmd, (char *)fmds);
  QDP_String *rmd = QDP_string_create();
  QDP_string_set(rmd, (char *)rmds);
  QDP_Writer *qw = QDP_open_write_L(g->qlat, fmd, (char*)fn, QDP_SINGLEFILE);
  int nd = g->nd;
  if(precision[0]=='F') {
#if QDP_Precision == 'F'
    QDP_vwrite_M(qw, rmd, g->links, nd);
#else
    QDP_F_ColorMatrix *tm[nd];
    for(int i=0; i<nd; i++) {
      tm[i] = QDP_F_create_M_L(g->qlat);
      QDP_FD_M_eq_M(tm[i], g->links[i], QDP_all_L(g->qlat));
    }
    QDP_F_vwrite_M(qw, rmd, tm, nd);
    for(int i=0; i<nd; i++) {
      QDP_F_destroy_M(tm[i]);
    }
#endif
  } else { // record precision is 'D'
#if QDP_Precision == 'F'
    QDP_D_ColorMatrix *tm[nd];
    for(int i=0; i<nd; i++) {
      tm[i] = QDP_D_create_M_L(g->qlat);
      QDP_DF_M_eq_M(tm[i], g->links[i], QDP_all_L(g->qlat));
    }
    QDP_D_vwrite_M(qw, rmd, tm, nd);
    for(int i=0; i<nd; i++) {
      QDP_D_destroy_M(tm[i]);
    }
#else
    QDP_vwrite_M(qw, rmd, g->links, nd);
#endif
  }
  QDP_close_write(qw);
  QDP_string_destroy(fmd);
  QDP_string_destroy(rmd);
  //dt += QDP_time();
  //printf0(" saved in %g seconds\n", dt);
  return 0;
#undef NC
}

// set from another field
static int
qopqdp_gauge_set(lua_State *L)
{
  BEGIN_ARGS;
  GET_GAUGE(g1);
  GET_GAUGE(g2);
  END_ARGS;
  if(g1->qlat==g2->qlat) {
    for(int i=0; i<g1->nd; i++) {
      QDP_M_eq_M(g1->links[i], g2->links[i], QDP_all_L(g1->qlat));
    }
  } else {
    // copy with truncation/replication
    QDP_Lattice *rlat = g1->qlat;
    QDP_Lattice *slat = g2->qlat;
    int rnd = QDP_ndim_L(rlat);
    int snd = QDP_ndim_L(slat);
    printf0("rls:");
    for(int i=0; i<rnd; i++) printf0(" %i", QDP_coord_size_L(rlat,i));
    printf0("\n");
    printf0("sls:");
    for(int i=0; i<snd; i++) printf0(" %i", QDP_coord_size_L(slat,i));
    printf0("\n");
    int roff[rnd], rlen[rnd], sdir[rnd], soff[snd];
    for(int i=0; i<rnd; i++) {
      roff[i] = 0;
      rlen[i] = QDP_coord_size_L(rlat,i);
      sdir[i] = i;
    }
    for(int i=0; i<snd; i++) {
      soff[i] = 0;
    }
    for(int s=0; ; s++) {
      printf0("s: %i\n", s);
      QDP_Subset *sub;
      QDP_Shift shift;
      qhmc_qopqdp_getCopyHyper(&shift, &sub, rlat, roff, rlen, sdir, slat, soff, s);
      if(sub==NULL) break;
      printf0("created subset and map\n");
      for(int mu=0; mu<rnd; mu++) {
	QDP_M_eq_sM(g1->links[mu], g2->links[mu], shift, QDP_forward, sub[0]);
      }
      printf0("finished shift\n");
      QDP_destroy_shift(shift);
      //printf0("destroyed map\n");
      QDP_destroy_subset(sub);
      //printf0("destroyed subset\n");
    }
  }
  return 0;
}

// create new copy
static int
qopqdp_gauge_copy(lua_State *L)
{
  //#define NC QDP_get_nc(g1->links[0])
  BEGIN_ARGS;
  GET_GAUGE(g1);
  OPT_INT(precision, QOP_PrecisionInt);
  END_ARGS;
  if(precision==QOP_PrecisionInt) {
    gauge_t *g2 = qopqdp_gauge_create(L, g1->nc, g1->lat);
    for(int i=0; i<g1->nd; i++) {
      QDP_M_eq_M(g2->links[i], g1->links[i], QDP_all_L(g1->qlat));
    }
  } else {
    gaugeO_t *g2 = qopqdp_gaugeO_create(L, g1->nc, g1->lat);
    for(int i=0; i<g1->nd; i++) {
      QDPOP(M_eq_M)(g2->links[i], g1->links[i], QDP_all_L(g1->qlat));
    }
  }
  return 1;
  //#undef NC
}

static void
checkU(NCPROT QLA_ColorMatrix(*m), int idx, void *args)
{
  double *devs = (double *) args;
  QLA_ColorMatrix(m2);
  QLA_M_eq_Ma_times_M(&m2, m, m);
  for(int i=0; i<QLA_Nc; i++) QLA_c_meq_r(QLA_elem_M(m2,i,i), 1);
  QLA_Real dev;
  QLA_r_eq_norm2_M(&dev, &m2);
  devs[0] += dev;
  if(dev>devs[1]) devs[1] = dev;
}

static int
qopqdp_gauge_checkU(lua_State *L)
{
#define NC QDP_get_nc(g->links[0])
  qassert(lua_gettop(L)==1);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  double devs[2];
  devs[0] = devs[1] = 0;
  for(int i=0; i<g->nd; i++) {
    QDP_M_eq_funcia(g->links[i], checkU, devs, QDP_all_L(g->qlat));
  }
  devs[0] /= g->nd*QDP_volume();
  devs[0] = sqrt(devs[0]/(2*QLA_Nc*QLA_Nc));
  devs[1] = sqrt(devs[1]/(2*QLA_Nc*QLA_Nc));
  lua_pushnumber(L, devs[0]);
  lua_pushnumber(L, devs[1]);
  return 2;
#undef NC
}

static void
checkSU(NCPROT QLA_ColorMatrix(*m), int idx, void *args)
{
  double *devs = (double *) args;
  QLA_ColorMatrix(m2);
  QLA_M_eq_Ma_times_M(&m2, m, m);
  for(int i=0; i<QLA_Nc; i++) QLA_c_meq_r(QLA_elem_M(m2,i,i), 1);
  QLA_Real dev;
  QLA_r_eq_norm2_M(&dev, &m2);
  QLA_Complex c;
  QLA_C_eq_det_M(&c, m);
  QLA_c_meq_r(c, 1);
  dev += QLA_norm2_c(c);
  devs[0] += dev;
  if(dev>devs[1]) devs[1] = dev;
}

static int
qopqdp_gauge_checkSU(lua_State *L)
{
#define NC QDP_get_nc(g->links[0])
  qassert(lua_gettop(L)==1);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  double devs[2];
  devs[0] = devs[1] = 0;
  for(int i=0; i<g->nd; i++) {
    QDP_M_eq_funcia(g->links[i], checkSU, devs, QDP_all_L(g->qlat));
  }
  devs[0] /= g->nd*QDP_volume();
  devs[0] = sqrt(devs[0]/(2*QLA_Nc*QLA_Nc+2));
  devs[1] = sqrt(devs[1]/(2*QLA_Nc*QLA_Nc+2));
  lua_pushnumber(L, devs[0]);
  lua_pushnumber(L, devs[1]);
  return 2;
#undef NC
}

void
qopqdp_makeSU(NCPROT QLA_ColorMatrix(*m), int idx, void *args)
{
  QLA_D_Complex d1,d2;
  QLA_ColorMatrix(m1);
  QLA_ColorMatrix(m2);
  QLA_Complex c1;
  QLA_M_eq_Ma_times_M(&m2, m, m);
  QLA_M_eq_sqrt_M(&m1, &m2);
  QLA_M_eq_inverse_M(&m2, &m1);
  QLA_M_eq_M_times_M(&m1, m, &m2);
  QLA_C_eq_det_M(&c1, &m1);
  QLA_c_eq_r_plus_ir(d1, QLA_real(c1), QLA_imag(c1));
  //d2 = QLA_clog(&d1);
  QLA_D_C_eq_clog_C(&d2, &d1);
  QLA_c_eq_r_times_c(d1, -1./QLA_Nc, d2);
  //d2 = QLA_cexp(&d1);
  QLA_D_C_eq_cexp_C(&d2, &d1);
  QLA_c_eq_r_plus_ir(c1, QLA_real(d2), QLA_imag(d2));
  QLA_M_eq_C_times_M(m, &c1, &m1);
}

static int
qopqdp_gauge_makeSU(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  for(int i=0; i<g->nd; i++) {
    QDP_M_eq_funcia(g->links[i], qopqdp_makeSU, NULL, QDP_all_L(g->qlat));
  }
  return 0;
}

static void
get_gauge_coeffs(lua_State *L, QOP_gauge_coeffs_t *coeffs, int idx)
{
  qassert(lua_type(L,idx)==LUA_TTABLE);
  *coeffs = QOP_GAUGE_COEFFS_ZERO;
  lua_pushnil(L);
  while(lua_next(L,idx)) {
    lua_pushvalue(L, -2);
    const char *s = lua_tostring(L, -1);
    if(strcmp(s,"plaq")==0) coeffs->plaquette = lua_tonumber(L,-2);
    if(strcmp(s,"rect")==0) coeffs->rectangle = lua_tonumber(L,-2);
    if(strcmp(s,"pgm")==0) coeffs->parallelogram = lua_tonumber(L,-2);
    if(strcmp(s,"adjplaq")==0) coeffs->adjoint_plaquette = lua_tonumber(L,-2);
    lua_pop(L, 2);
  }
}

// 1: gauge
// 2: coeffs
// 3: (optional) xi0
static int
qopqdp_gauge_action(lua_State *L)
{
  int nargs = lua_gettop(L);
  qassert(nargs>=2&&nargs<=3);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  QOP_gauge_coeffs_t coeffs;
  get_gauge_coeffs(L, &coeffs, 2);
  QLA_Real ixi0=1, xi0=luaL_optnumber(L, 3, 1);
  if(xi0!=1) {
    ixi0 = 1/xi0;
    QDP_M_eq_r_times_M(g->links[g->nd-1], &xi0, g->links[g->nd-1], QDP_all_L(g->qlat));
  }
  QOP_info_t info;
  QLA_Real acts, actt;
  QOP_GaugeField *qg = QOP_create_G_from_qdp(g->links);
  QOP_symanzik_1loop_gauge_action(&info, qg, &acts, &actt, &coeffs);
  QOP_destroy_G(qg);
  if(xi0!=1) {
    QDP_M_eq_r_times_M(g->links[g->nd-1], &ixi0, g->links[g->nd-1], QDP_all_L(g->qlat));
  }
  double ma = 0;
  //ma = get_measure_action(g);
  lua_pushnumber(L, ixi0*acts);
  lua_pushnumber(L, ixi0*actt);
  lua_pushnumber(L, ma);
  return 3;
}

// 1: gauge
// 2: force
// 3: coeffs
// 4: (optional) beta
// 5: (optional) xi0
static int
qopqdp_gauge_force(lua_State *L)
{
  int nargs = lua_gettop(L);
  qassert(nargs>=3&&nargs<=5);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  force_t *f = qopqdp_force_check(L, 2);
  for(int i=0; i<f->nd; i++) QDP_M_eq_zero(f->force[i], QDP_all_L(g->qlat));
  QOP_gauge_coeffs_t coeffs;
  get_gauge_coeffs(L, &coeffs, 3);
  double beta = luaL_optnumber(L, 4, 1);
  QLA_Real ixi0=1, xi0=luaL_optnumber(L, 5, 1);
  if(xi0!=1) {
    ixi0 = 1/xi0;
    beta *= ixi0;
    QDP_M_eq_r_times_M(g->links[g->nd-1], &xi0, g->links[g->nd-1], QDP_all_L(g->qlat));
  }
  QOP_info_t info;
  QOP_GaugeField *qg = QOP_create_G_from_qdp(g->links);
  QOP_Force *qf = QOP_create_F_from_qdp(f->force);
  switch(g->nc) {
#ifdef HAVE_NC3
  case 3:
    QOP_3_symanzik_1loop_gauge_force(&info, (QOP_3_GaugeField*)qg,
				     (QOP_3_Force*)qf, &coeffs, beta);
    break;
#endif
  default: 
   QOP_symanzik_1loop_gauge_force(&info, qg, qf, &coeffs, beta);
  }
  QOP_extract_F_to_qdp(f->force, qf);
  QOP_destroy_G(qg);
  QOP_destroy_F(qf);
  if(xi0!=1) {
    QDP_M_eq_r_times_M(g->links[g->nd-1], &ixi0, g->links[g->nd-1], QDP_all_L(g->qlat));
  }
  f->time = info.final_sec;
  f->flops = info.final_flop;
  return 0;
}

static int
qopqdp_gauge_update(lua_State *L)
{
#define NC QDP_get_nc(g->links[0])
  qassert(lua_gettop(L)==3);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  force_t *f = qopqdp_force_check(L, 2);
  int nd = f->nd;
  double eps[nd];
  if(lua_type(L,3)==LUA_TTABLE) {
    qhmc_get_double_array(L, 3, nd, eps);
  } else {
    double s = luaL_checknumber(L, 3);
    for(int i=0; i<nd; i++) eps[i] = s;
  }
#if 0
  {
    QLA_ColorMatrix *cm = QDP_expose_M(f->force[0]);
    printf("%g\t%g\n", eps[0], QLA_imag(QLA_elem_M(*cm,0,0)));
    QDP_reset_M(f->force[0]);
  }
#endif
  //get_gauge_links(g);
  QDP_ColorMatrix *m1 = QDP_create_M_L(g->qlat);
  QDP_ColorMatrix *m2 = QDP_create_M_L(g->qlat);
  for(int i=0; i<nd; i++) {
    QLA_Real teps = eps[i];
    QDP_M_eq_r_times_M(m1, &teps, f->force[i], QDP_all_L(g->qlat));
    //QDP_M_eq_exp_M(m2, m1, QDP_all_L(m2));
    QDP_M_eq_expTA_M(m2, m1, QDP_all_L(g->qlat));  // assume traceless antihermitian
    QDP_M_eq_M_times_M(m1, m2, g->links[i], QDP_all_L(g->qlat));
    QDP_M_eq_M(g->links[i], m1, QDP_all_L(g->qlat));
#ifdef QHMC_REPRO_UNIFORM
    int err;
#define chk(x, i, s) \
    err = check_uniform_M((x)[i], s); \
    if(err) { \
      fprintf(stderr, "%i repro error %s %i index %i\n", QDP_this_node, #x, i, err-1); \
      QDP_abort(1); \
    }
    chk(f->force, i, QDP_even);
    chk(f->force, i, QDP_odd);
    chk(g->links, i, QDP_even);
    chk(g->links, i, QDP_odd);
#endif
  }
  QDP_destroy_M(m1);
  QDP_destroy_M(m2);
  return 0;
#undef NC
}

// 1: gauge
// 2: nrep (# of repetitions)
// 3: nhb (# heatbath updates per repetition)
// 4: nor (# overrelaxation updates per repetition)
// 5: beta
// 6: coeffs
// 7: xi0 (unused)
static int
qopqdp_gauge_heatbath(lua_State *L)
{
  BEGIN_ARGS;
  GET_GAUGE(g);
  GET_INT(nrep);
  GET_INT(nhb);
  GET_INT(nor);
  GET_DOUBLE(beta);
  GET_GAUGE_COEFFS(coeffs);
  //OPT_DOUBLE(xi0, 1);
  END_ARGS;
  QOP_info_t info;
  QOP_symanzik_1loop_gauge_heatbath_qdp(&info, g->links, beta, &coeffs,
					qopqdp_srs, nrep, nhb, nor);
  return 0;
}

// calculate generic loop
// 1: gauge field
// 2: path (+(1+mu): from forward, -(1+mu): from backward)
// 3: (opt) subset/subsets
static int
qopqdp_gauge_loop(lua_State *L)
{
#define NC QDP_get_nc(g->links[0])
  BEGIN_ARGS;
  GET_GAUGE(g);
  GET_TABLE_LEN_INDEX(plen,pidx);
  OPT_SUBSETS(subs, ns, g->lat, QDP_all_and_empty_L(g->qlat), 1);
  END_ARGS;
  int dirs[plen]; qhmc_get_int_array(L, pidx, plen, dirs);
  QDP_Lattice *qlat = g->qlat;
  QDP_Subset all = QDP_all_L(qlat);
  QDP_Shift *neighbor = QDP_neighbor_L(qlat);

  QLA_Complex z;
  QDP_ColorMatrix *m[2], *temp[2];
  m[0] = QDP_create_M_L(qlat);
  m[1] = QDP_create_M_L(qlat);
  temp[0] = QDP_create_M_L(qlat);
  temp[1] = QDP_create_M_L(qlat);

  int k = 0;
  int intemp = 0;
  QLA_c_eq_r(z, 1);
  QDP_M_eq_c(m[k], &z, all);
  for(int i=0; i<plen; i++) {
    int d = dirs[i];
    int mu = abs(d) - 1;
    //printf0("i: %i  d: %i  mu: %i\n", i, d, mu);
    if(mu<0 || mu>=g->nd) {
      printf0("%s invalid index[%i] %i (nd=%i)\n", __func__, i, d, g->nd);
      lua_pushnil(L);
      return 1;
    }
    if(d<0) { // shift from backward
      if(intemp) {
	QDP_M_eq_Ma_times_M(m[1-k], g->links[mu], temp[k], all);
	//QDP_discard_M(temp[k]);
      } else {
	QDP_M_eq_Ma_times_M(m[1-k], g->links[mu], m[k], all);
      }
      k = 1-k;
      QDP_M_eq_sM(temp[k], m[k], neighbor[mu], QDP_backward, all);
      intemp = 1;
    } else { // shift from forward
      if(intemp) {
	QDP_M_eq_M(m[1-k], temp[k], all);
	//QDP_discard_M(temp[k]);
	k = 1-k;
      }
      //printf("k: %i\n", k);
      //TRACE;
      QDP_M_eq_sM(temp[k], m[k], neighbor[mu], QDP_forward, all);
      //TRACE;
      QDP_M_eq_M_times_M(m[1-k], g->links[mu], temp[k], all);
      //TRACE;
      k = 1-k;
      intemp = 0;
    }
  }
  if(intemp) {
    QDP_M_eq_M(m[1-k], temp[k], all);
    QDP_discard_M(temp[k]);
    k = 1-k;
  }

  QLA_Complex qc[ns];
  QDP_Complex *c = QDP_create_C_L(qlat);
  QDP_C_eq_trace_M(c, m[k], all);
  QDP_c_eq_sum_C_multi(qc, c, subs, ns);
  QDP_destroy_C(c);
  double f = 1.0/(QLA_Nc*QDP_volume_L(qlat));
  qhmc_complex_t cc[ns];
  for(int i=0; i<ns; i++) {
    cc[i].r = f*QLA_real(qc[i]);
    cc[i].i = f*QLA_imag(qc[i]);
  }
  if(ns==1) {
    qhmc_complex_create(L, cc[0].r, cc[0].i);
  } else {
    qhmc_push_complex_array(L, ns, cc);
  }

  QDP_destroy_M(m[0]);
  QDP_destroy_M(m[1]);
  QDP_destroy_M(temp[0]);
  QDP_destroy_M(temp[1]);

  return 1;
#undef NC
}

// calculate s4 broken phase observables
// Based on "meas_plaq" in milc
static int
qopqdp_gauge_s4_gauge_observables(lua_State *L)
{
#define NC QDP_get_nc(g->links[0])
  // All we expect is the gauge field.
  qassert(lua_gettop(L)==1);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  lattice_t *lat = g->lat;

  // We need to construct the plaquette!
  QDP_ColorMatrix *plaq, *temp, *temp2;
  QDP_Real* retrplaq;
	
  // For creating subsections.
  QDP_Subset *dirsubset;
	
  // According to "meas_plaq," we need these too.
  // 0 = spacetime  (0 == even, 1 == odd)
  // 1 = contains x
  // 2 = contains y
  // 3 = contains z
  // 4 = all coords odd, all coords even (opposite convention)
  QLA_Real plaq_ss = 0, plaq_e[5], plaq_o[5];
  QLA_Real interm[3];
  int i, dir1, dir2;
	
  plaq = QDP_create_M();
  temp = QDP_create_M();
  temp2 = QDP_create_M();
  retrplaq = QDP_create_R();
	
  for (i = 0; i < 5; i++)
    {
      plaq_e[i] = plaq_o[i] = 0;
    }
	
  // Get that plaquette! Hardcoded for 4D.
  for (dir1 = 1; dir1 < 4; dir1++)
    {
      for (dir2 = 0; dir2 < dir1; dir2++)
	{
	  // Get dat plaq!
	  QDP_M_eq_sM(temp, g->links[dir1], QDP_neighbor[dir2], QDP_forward, QDP_all);
	  
	  QDP_M_eq_M_times_M(plaq, g->links[dir2], temp, QDP_all);
	  
	  QDP_M_eq_sM(temp, g->links[dir2], QDP_neighbor[dir1], QDP_forward, QDP_all);
	  
	  QDP_M_eq_M_times_Ma(temp2, plaq, temp, QDP_all);
	  
	  QDP_M_eq_M_times_Ma(plaq, temp2, g->links[dir1], QDP_all);
	  
	  QDP_R_eq_re_trace_M(retrplaq, plaq, QDP_all);
			
	  // We have the plaq!
	  
	  if (dir1 == 3) // time-space
	    {
	      dirsubset = qhmcqdp_get_eodir(lat, 3);
	      QDP_r_eq_sum_R_multi(interm, retrplaq, dirsubset, 2);
	      // interm[0] is the sum over values where t%2==0
	      // interm[1] is the sum over values where t%2==1
	      plaq_e[0] += interm[0];
	      plaq_o[0] += interm[1];
	    }
	  else // spatial
	    {
	      QDP_r_eq_sum_R(interm, retrplaq, QDP_all);
	      plaq_ss += interm[0];
	    }
	  
	  if (dir1 == 0 || dir2 == 0) // x dep.
	    {
	      dirsubset = qhmcqdp_get_eodir(lat, 0);
	      QDP_r_eq_sum_R_multi(interm, retrplaq, dirsubset, 2);
			    
	      // interm[0] is the sum over values where x%2==0
	      // interm[1] is the sum over values where x%2==1
	      plaq_e[1] += interm[0];
	      plaq_o[1] += interm[1];
	    }
	  
	  if (dir1 == 1 || dir2 == 1) // y dep.
	    {
	      dirsubset = qhmcqdp_get_eodir(lat, 1);
	      QDP_r_eq_sum_R_multi(interm, retrplaq, dirsubset, 2);
	      
	      // interm[0] is the sum over values where y%2==0
	      // interm[1] is the sum over values where y%2==1
	      plaq_e[2] += interm[0];
	      plaq_o[2] += interm[1];
	    }
	  
	  if (dir1 == 2 || dir2 == 2) // z dep.
	    {
	      dirsubset = qhmcqdp_get_eodir(lat, 2);
	      QDP_r_eq_sum_R_multi(interm, retrplaq, dirsubset, 2);
	      
	      // interm[0] is the sum over values where z%2==0
	      // interm[1] is the sum over values where z%2==1
	      plaq_e[3] += interm[0];
	      plaq_o[3] += interm[1];
	    }
	  
	  // The last weird case. All coordinates even or all
	  // coordinates odd.
	  dirsubset = qhmcqdp_get_eodir(lat, 4);
	  
	  QDP_r_eq_sum_R_multi(interm, retrplaq, dirsubset, 3);
	  
	  plaq_e[4] += interm[1]; // Yes, I know it's switched.
	  plaq_o[4] += interm[0]; // This matches MILC.
	}
    }
  
  // Normalize.
  plaq_ss /= ((double)3*3*QDP_volume());
  for (i = 0; i < 4; i++)
    {
      plaq_e[i] /= ((double)3*3*QDP_volume())/2;
		plaq_o[i] /= ((double)3*3*QDP_volume())/2;
    }
  plaq_e[4] /= ((double)3*3*QDP_volume())/8;
  plaq_o[4] /= ((double)3*3*QDP_volume())/8;
  
  // Clean up!
  QDP_destroy_M(plaq);
  QDP_destroy_M(temp);
  QDP_destroy_M(temp2);
  QDP_destroy_R(retrplaq);
  
  // Return things.
  lua_pushnumber(L, plaq_ss);
  qhmc_push_real_array(L, 5, plaq_e);
  qhmc_push_real_array(L, 5, plaq_o);

  return 3;
#undef NC
}

static struct luaL_Reg gauge_reg[] = {
  { "__gc",     qopqdp_gauge_gc },
  { "__len",    qopqdp_gauge_len },
  { "__call",   qopqdp_gauge_call },
  { "nc",       qopqdp_gauge_nc },
  { "lattice",  qopqdp_gauge_lattice },
  { "unit",     qopqdp_gauge_unit },
  { "random",   qopqdp_gauge_random },
  { "load",     qopqdp_gauge_load },
  { "save",     qopqdp_gauge_save },
  { "set",      qopqdp_gauge_set },
  { "copy",     qopqdp_gauge_copy },
  { "checkU",   qopqdp_gauge_checkU },
  { "checkSU",  qopqdp_gauge_checkSU },
  { "makeSU",   qopqdp_gauge_makeSU },
  { "action",   qopqdp_gauge_action },
  { "force",    qopqdp_gauge_force },
  { "update",   qopqdp_gauge_update },
  { "heatbath", qopqdp_gauge_heatbath },
  { "loop",     qopqdp_gauge_loop },
  { "coulomb",  qopqdp_gauge_coulomb },
  { "s4Gauge",  qopqdp_gauge_s4_gauge_observables }, // ESW addition 12/18/2013
  { NULL, NULL}
};

// different fwd & bck links

gauge_t *
qopqdp_gauge_create(lua_State *L, int nc, lattice_t *lat)
{
#define NC nc
  if(lat==NULL) lat = qopqdp_get_default_lattice(L);
  if(nc==0) nc = lat->defaultNc;
  int nd = QDP_ndim_L(lat->qlat);
  gauge_t *g = lua_newuserdata(L,sizeof(gauge_t)+nd*sizeof(QDP_ColorMatrix*));
  g->nd = nd;
  g->lat = lat;
  g->qlat = lat->qlat;
  g->nc = nc;
  for(int i=0; i<nd; i++) {
    g->links[i] = QDP_create_M_L(lat->qlat);
  }
  if(luaL_newmetatable(L, gmtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, gauge_reg);
  }
  lua_setmetatable(L, -2);
  return g;
#undef NC
}
