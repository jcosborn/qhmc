#include "qhmc_qopqdp_common.h"
#include <string.h>
#include <math.h>

static char *gmtname = "qopqdp.gauge";

gauge_t *
qopqdp_gauge_check(lua_State *L, int idx)
{
  luaL_checkudata(L, idx, gmtname);
  gauge_t *g = lua_touserdata(L, idx);
#if 0
  int hasmt = lua_getmetatable(L, idx);
  qassert(hasmt==1);
  luaL_getmetatable(L, gmtname);
  int eq = lua_equal(L, -1, -2);
  qassert(eq==1);
  lua_pop(L, 2);
#endif
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
    if(g->lie) QDP_destroy_M(g->lie[i]);
  }
  free(g->links);
  if(g->lie) free(g->lie);
}

static int
qopqdp_gauge_gc(lua_State *L)
{
  qopqdp_gauge_free(L, -1);
  return 0;
}

static void
get_gauge_links(gauge_t *g)
{
  QDP_ColorMatrix *m = QDP_create_M();
  for(int i=0; i<g->nd; i++) {
    QDP_M_eq_i_M(m, g->lie[i], QDP_all);
    QDP_M_eq_exp_M(g->links[i], m, QDP_all);
  }
  QDP_destroy_M(m);
}

static void
get_gauge_lie(gauge_t *g)
{
  QDP_ColorMatrix *m = QDP_create_M();
  for(int i=0; i<g->nd; i++) {
    QDP_M_eq_log_M(m, g->links[i], QDP_all);
    QDP_M_eqm_i_M(g->lie[i], m, QDP_all);
  }
  QDP_destroy_M(m);
}

static int
qopqdp_gauge_unit(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  gauge_t *g = qopqdp_gauge_check(L, -1);
  QLA_Complex z;
  QLA_c_eq_r(z, 1);
  for(int i=0; i<g->nd; i++) {
    if(g->lie) QDP_M_eq_zero(g->lie[i], QDP_all);
    QDP_M_eq_c(g->links[i], &z, QDP_all);
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
  QLA_M_eq_sqrt_M(&m3, &m2);
  QLA_M_eq_inverse_M(&m2, &m3);
  QLA_M_eq_M_times_M(&m3, &m1, &m2);
  QLA_C_eq_det_M(&c, &m3);
  QLA_DF_c_eq_c(d1, c);
  d2 = QLA_clog(&d1);
  QLA_c_eq_r_times_c(d1, -1./QLA_Nc, d2);
  d2 = QLA_cexp(&d1);
  QLA_FD_c_eq_c(c, d2);
  QLA_M_eq_C_times_M(m, &c, &m3);
}

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

static int
qopqdp_gauge_random(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  if(g->lie) {
    for(int i=0; i<g->nd; i++) {
      QDP_M_eq_gaussian_S(g->lie[i], qopqdp_srs, QDP_all);
      QDP_M_eq_funcia(g->lie[i], make_herm, NULL, QDP_all);
    }
    get_gauge_links(g);
  } else {
    QLA_RandomState *qrs = QDP_expose_S(qopqdp_srs);
    for(int i=0; i<g->nd; i++) {
      QDP_M_eq_funcia(g->links[i], gauge_random, qrs, QDP_all);
    }
    QDP_reset_S(qopqdp_srs);
  }
  return 0;
}

/* get QIO record precision */
static int
get_prec(QDP_Reader *qr)
{
  QIO_RecordInfo *ri = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);
  QDP_String *md = QDP_string_create();
  QDP_read_qio_record_info(qr, ri, md);
  int prec = *QIO_get_precision(ri);
  QIO_destroy_record_info(ri);
  QDP_string_destroy(md);
  return prec;
}

static int
qopqdp_gauge_load(lua_State *L)
{
  qassert(lua_gettop(L)==2);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  const char *fn = luaL_checkstring(L, 2);

  QDP_set_read_group_size(8);
  printf0("loading lattice file %s\n", fn);
  double dt = -QDP_time();
  QDP_String *md = QDP_string_create();
  QDP_Reader *qr = QDP_open_read(md, (char*)fn);
  int nd = QDP_ndim();
  if(get_prec(qr)=='F') {
#if QDP_Precision == 'F'
    QDP_F_vread_M(qr, md, g->links, nd);
#else
    QDP_F_ColorMatrix *tm[nd];
    for(int i=0; i<nd; i++) {
      tm[i] = QDP_F_create_M();
    }
    QDP_F_vread_M(qr, md, tm, nd);
    for(int i=0; i<nd; i++) {
      QDP_DF_M_eq_M(g->links[i], tm[i], QDP_all);
      QDP_F_destroy_M(tm[i]);
    }
#endif
  } else { // record precision is 'D'
#if QDP_Precision == 'F'
    QDP_D_ColorMatrix *tm[nd];
    for(int i=0; i<nd; i++) {
      tm[i] = QDP_D_create_M();
    }
    QDP_D_vread_M(qr, md, tm, nd);
    for(int i=0; i<nd; i++) {
      QDP_FD_M_eq_M(g->links[i], tm[i], QDP_all);
      QDP_D_destroy_M(tm[i]);
    }
#else
    QDP_D_vread_M(qr, md, g->links, nd);
#endif
  }
  QDP_close_read(qr);
  QDP_string_destroy(md);
  if(g->lie) get_gauge_lie(g);
  dt += QDP_time();
  printf0(" loaded in %g seconds\n", dt);
  return 0;
}

static int
qopqdp_gauge_save(lua_State *L)
{
  int nargs = lua_gettop(L);
  qassert(nargs==3 || nargs==4);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  const char *fn = luaL_checkstring(L, 2);
  const char *mds = luaL_checkstring(L, 3);
  char prec = 'F';
  if(nargs>3) {
    prec = *luaL_checkstring(L, 4);
  }

  QDP_set_write_group_size(8);
  printf0("saving lattice file %s\n", fn);
  double dt = -QDP_time();
  QDP_String *md = QDP_string_create();
  QDP_string_set(md, (char *)mds);
  QDP_Writer *qw = QDP_open_write(md, (char*)fn, QDP_SINGLEFILE);
  int nd = QDP_ndim();
  if(prec=='F') {
#if QDP_Precision == 'F'
    QDP_vwrite_M(qw, md, g->links, nd);
#else
    QDP_F_ColorMatrix *tm[nd];
    for(int i=0; i<nd; i++) {
      tm[i] = QDP_F_create_M();
      QDP_FD_M_eq_M(tm[i], g->links[i], QDP_all);
    }
    QDP_F_vwrite_M(qw, md, tm, nd);
    for(int i=0; i<nd; i++) {
      QDP_F_destroy_M(tm[i]);
    }
#endif
  } else { // record precision is 'D'
#if QDP_Precision == 'F'
    QDP_D_ColorMatrix *tm[nd];
    for(int i=0; i<nd; i++) {
      tm[i] = QDP_D_create_M();
      QDP_D_M_eq_M(tm[i], g->links[i], QDP_all);
    }
    QDP_D_vwrite_M(qw, md, tm, nd);
    for(int i=0; i<nd; i++) {
      QDP_D_destroy_M(tm[i]);
    }
#else
    QDP_vwrite_M(qw, md, g->links, nd);
#endif
  }
  QDP_close_write(qw);
  QDP_string_destroy(md);
  dt += QDP_time();
  printf0(" saved in %g seconds\n", dt);
  return 0;
}

static int
qopqdp_gauge_set(lua_State *L)
{
  int nargs = lua_gettop(L);
  qassert(nargs==2);
  gauge_t *g1 = qopqdp_gauge_check(L, 1);
  gauge_t *g2 = qopqdp_gauge_check(L, 2);
  qassert(g1->nd==g2->nd);
  for(int i=0; i<g1->nd; i++) {
    QDP_M_eq_M(g1->links[i], g2->links[i], QDP_all);
    if(g1->lie&&g2->lie) QDP_M_eq_M(g1->lie[i], g2->lie[i], QDP_all);
  }
  return 0;
}

static int
qopqdp_gauge_copy(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  gauge_t *g1 = qopqdp_gauge_check(L, 1);
  gauge_t *g2 = qopqdp_gauge_create(L);
  g2->nd = g1->nd;
  for(int i=0; i<g2->nd; i++) {
    QDP_M_eq_M(g2->links[i], g1->links[i], QDP_all);
    if(g1->lie&&g2->lie) QDP_M_eq_M(g2->lie[i], g1->lie[i], QDP_all);
  }
  return 1;
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
  qassert(lua_gettop(L)==1);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  double devs[2];
  devs[0] = devs[1] = 0;
  for(int i=0; i<g->nd; i++) {
    QDP_M_eq_funcia(g->links[i], checkU, devs, QDP_all);
  }
  devs[0] /= g->nd*QDP_volume();
  devs[0] = sqrt(devs[0]/(2*QLA_Nc*QLA_Nc));
  devs[1] = sqrt(devs[1]/(2*QLA_Nc*QLA_Nc));
  lua_pushnumber(L, devs[0]);
  lua_pushnumber(L, devs[1]);
  return 2;
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
  qassert(lua_gettop(L)==1);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  double devs[2];
  devs[0] = devs[1] = 0;
  for(int i=0; i<g->nd; i++) {
    QDP_M_eq_funcia(g->links[i], checkSU, devs, QDP_all);
  }
  devs[0] /= g->nd*QDP_volume();
  devs[0] = sqrt(devs[0]/(2*QLA_Nc*QLA_Nc+2));
  devs[1] = sqrt(devs[1]/(2*QLA_Nc*QLA_Nc+2));
  lua_pushnumber(L, devs[0]);
  lua_pushnumber(L, devs[1]);
  return 2;
}

static void
makeSU(NCPROT QLA_ColorMatrix(*m), int idx, void *args)
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
  QLA_C_eq_clog_C(&d2, &d1);
  QLA_c_eq_r_times_c(d1, -1./QLA_Nc, d2);
  //d2 = QLA_cexp(&d1);
  QLA_C_eq_cexp_C(&d2, &d1);
  QLA_c_eq_r_plus_ir(c1, QLA_real(d2), QLA_imag(d2));
  QLA_M_eq_C_times_M(m, &c1, &m1);
}

static int
qopqdp_gauge_makeSU(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  if(g->lie) {
#if 0
    {
      printf0("initial mat:\n");
      QLA_ColorMatrix *cm = QDP_expose_M(g->lie[0]);
      for(int i=0; i<QLA_Nc; i++) {
	for(int j=0; j<QLA_Nc; j++) {
	  QLA_Complex *z = &QLA_elem_M(*cm, i, j);
	  printf0(" lie%i%i\t%g\t%g\n", i, j, QLA_real(*z), QLA_imag(*z));
	}
      }
      QDP_reset_M(g->lie[0]);
    }
#endif
    for(int i=0; i<g->nd; i++) {
      QDP_M_eq_funcia(g->lie[i], make_herm, NULL, QDP_all);
    }
#if 0
    {
      QLA_Real max = 0;
      for(int i=0; i<g->nd; i++) {
	QLA_Real tmax = infnorm_M(g->lie[i], QDP_all);
	if(tmax>max) max = tmax;
      }
      printf0("max = %g\n", max);
    }
    get_gauge_links(g);
    get_gauge_lie(g);
    for(int i=0; i<g->nd; i++) {
      QDP_M_eq_funcia(g->lie[i], make_herm, NULL, QDP_all);
    }
    {
      QLA_Real max = 0;
      for(int i=0; i<g->nd; i++) {
	QLA_Real tmax = infnorm_M(g->lie[i], QDP_all);
	if(tmax>max) max = tmax;
      }
      printf0("max = %g\n", max);
    }
#endif
    get_gauge_links(g);
  } else {
    for(int i=0; i<g->nd; i++) {
      QDP_M_eq_funcia(g->links[i], makeSU, NULL, QDP_all);
    }
  }
  return 0;
}

#if 0
static double
getfieldnum(lua_State *L, int idx, char *s, double def)
{
  lua_getfield(L, idx, s);
  if(!lua_isnil(L,-1)) {
    qassert(lua_isnumber(L,-1));
    def = lua_tonumber(L, -1);
  }
  lua_pop(L, 1);
  return def;
}
#endif

static void
get_gauge_coeffs(lua_State *L, QOP_gauge_coeffs_t *coeffs, int idx)
{
  qassert(lua_type(L,idx)==LUA_TTABLE);
  coeffs->plaquette = 0;
  coeffs->rectangle = 0;
  coeffs->parallelogram = 0;
  lua_pushnil(L);
  while(lua_next(L,idx)) {
    lua_pushvalue(L, -2);
    const char *s = lua_tostring(L, -1);
    if(strcmp(s,"plaq")==0) coeffs->plaquette = lua_tonumber(L,-2);
    if(strcmp(s,"rect")==0) coeffs->rectangle = lua_tonumber(L,-2);
    if(strcmp(s,"pgm")==0) coeffs->parallelogram = lua_tonumber(L,-2);
    lua_pop(L, 2);
  }
}

#if 0
static double
sinc(double x)
{
  double y = 1;
  if(x!=0) y = sin(x)/x;
  return y;
}

static double
measure_action(QLA_ColorMatrix *m)
{
  // assume m is traceless Hermitian
  QLA_ColorMatrix m2;
  double p2, p3;
  QLA_R_eq_re_M_dot_M(&p2, m, m);
  QLA_M_eq_M_times_M(&m2, m, m);
  QLA_R_eq_re_M_dot_M(&p3, m, &m2);
  double a = 4.5*p3;
  double b = 1.5*p2;
  double c = sqrt(fabs(b*b*b-a*a));
  double r = sqrt(b);
  double t = atan2(c,a)/3.;
  double st = sin(t);
  double ct = sqrt(1-st*st);
  double y = r*st/sqrt(3.);
  double x = r*ct;
  //double l0 = 2*x/3;
  //double l1 = y-x/3;
  //double l2 = -y-x/3;
  double d0 = y;
  double d1 = 0.5*(x+y);
  double d2 = 0.5*(x-y);
  double z = sinc(d0)*sinc(d1)*sinc(d2);
  return -log(z*z);
  //return 0;
  //return p2;
}

static void
Qmeasure_action(QLA_ColorMatrix *m, int idx, void *args)
{
  double a = measure_action(m);
  *(double*)args += a;
}  

static double
get_measure_action(gauge_t *g)
{
  double a = 0;
  for(int i=0; i<g->nd; i++) {
    QDP_M_eq_funcia(g->lie[i], Qmeasure_action, &a, QDP_all);
  }
  return a;
}
#endif

static int
qopqdp_gauge_action(lua_State *L)
{
  qassert(lua_gettop(L)==2);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  QOP_gauge_coeffs_t coeffs;
  get_gauge_coeffs(L, &coeffs, 2);
  QOP_info_t info;
  QLA_Real acts, actt;
  QOP_GaugeField *qg = QOP_create_G_from_qdp(g->links);
  QOP_symanzik_1loop_gauge_action(&info, qg, &acts, &actt, &coeffs);
  QOP_destroy_G(qg);
  double ma = 0;
  //ma = get_measure_action(g);
  lua_pushnumber(L, acts);
  lua_pushnumber(L, actt);
  lua_pushnumber(L, ma);
  return 3;
}

#if 0
static void
update_link(QLA_ColorMatrix *g, QLA_ColorMatrix *g0, QLA_ColorMatrix *dg, double eps)
{
  //QLA_Complex ieps;
  //QLA_c_eq_r_plus_ir(ieps, 0, -eps);
  //QLA_M_eq_c_times_M_plus_M(g, &ieps, dg, g0);
  QLA_Real qeps = eps;
  QLA_M_eq_r_times_M_plus_M(g, &qeps, dg, g0);
}

static double
group_act(QLA_ColorMatrix *g, void *args)
{
  return measure_action(g);
}

static void
lie_force2(QLA_ColorMatrix *m, int idx, void *args)
{
  QLA_ColorMatrix *g0 = ((QLA_ColorMatrix *)args) + idx;
  QLA_ColorMatrix g, e0, e, im, f;
  double eps = 1e-6;
#define eim(y,x) QLA_M_eq_i_M(&im, x); QLA_M_eq_exp_M(y, &im)
  QLA_M_eq_zero(&f);
  eim(&e0, g0);
  for(int i=0; i<QLA_Nc; i++) {
    for(int j=i+1; j<QLA_Nc; j++) {
      QLA_Real zr, zi;

      QLA_M_eq_M(&g, g0);
      QLA_c_peq_r(QLA_elem_M(g,i,j), eps);
      QLA_c_peq_r(QLA_elem_M(g,j,i), eps);
      eim(&e, &g);
      QLA_M_meq_M(&e, &e0);
      QLA_R_eq_re_M_dot_M(&zr, &e, m);
      zr /= 2*eps;

      QLA_M_eq_M(&g, g0);
      QLA_c_peq_r_plus_ir(QLA_elem_M(g,i,j), 0, eps);
      QLA_c_peq_r_plus_ir(QLA_elem_M(g,j,i), 0, -eps);
      eim(&e, &g);
      QLA_M_meq_M(&e, &e0);
      QLA_R_eq_re_M_dot_M(&zi, &e, m);
      zi /= 2*eps;

      QLA_c_eq_r_plus_ir(QLA_elem_M(f,i,j), zr, zi);
      QLA_c_eq_r_plus_ir(QLA_elem_M(f,j,i), zr, -zi);
#if 0
      if(idx==0) {
	printf0("%i %i\t%g\t%g\n", i, j, zr, zi);
      }
#endif
    }
  }
  {
    QLA_Real zr;
    QLA_M_eq_M(&g, g0);
    QLA_c_peq_r(QLA_elem_M(g,0,0), eps);
    QLA_c_peq_r(QLA_elem_M(g,1,1), -eps);
    eim(&e, &g);
    QLA_M_meq_M(&e, &e0);
    QLA_R_eq_re_M_dot_M(&zr, &e, m);
    zr /= 2*eps;
    QLA_c_peq_r(QLA_elem_M(f,0,0), zr);
    QLA_c_peq_r(QLA_elem_M(f,1,1), -zr);
    //if(idx==0) printf("f3 = %g\n", zr);
  }
  {
    double s = sqrt(1./6.);
    QLA_Real zr;
    QLA_M_eq_M(&g, g0);
    QLA_c_peq_r(QLA_elem_M(g,0,0), s*eps);
    QLA_c_peq_r(QLA_elem_M(g,1,1), s*eps);
    QLA_c_peq_r(QLA_elem_M(g,2,2), -2*s*eps);
    eim(&e, &g);
    QLA_M_meq_M(&e, &e0);
    QLA_R_eq_re_M_dot_M(&zr, &e, m);
    zr *= s/eps;
    QLA_c_peq_r(QLA_elem_M(f,0,0), zr);
    QLA_c_peq_r(QLA_elem_M(f,1,1), zr);
    QLA_c_peq_r(QLA_elem_M(f,2,2), -2*zr);
  }

  QLA_M_eqm_i_M(m, &f);
  //QLA_Real s = 1;
  //QLA_M_eq_r_times_M(m, &s, m);
}

void
get_lie_force(force_t *f, gauge_t *g)
{
  for(int i=0; i<f->nd; i++) {
    QLA_ColorMatrix *m = QDP_expose_M(g->lie[i]);
    QDP_M_eq_funcia(f->force[i], lie_force2, m, QDP_all);
    QDP_reset_M(g->lie[i]);
  }
}

static void
measure_force(QLA_ColorMatrix *m, int idx, void *args)
{
  QLA_ColorMatrix f;
  QLA_ColorMatrix *g = ((QLA_ColorMatrix *)args) + idx;
  get_local_force(&f, g, update_link, group_act, NULL, 1e-6);
  QLA_M_peq_M(m, &f);
}

static void
get_measure_force(force_t *f, gauge_t *g)
{
  for(int i=0; i<f->nd; i++) {
    QLA_ColorMatrix *m = QDP_expose_M(g->lie[i]);
    QDP_M_eq_funcia(f->force[i], measure_force, m, QDP_all);
    QDP_reset_M(g->lie[i]);
  }
}
#endif

#if 0
static double
action(gauge_t *g, void *args)
{
  void **aa = (void **) args;
  QOP_gauge_coeffs_t *coeffs = (QOP_gauge_coeffs_t *) aa[0];
  double beta = * (double *) aa[1];
  QOP_info_t info;
  QLA_Real acts, actt;
  QOP_GaugeField *qg = QOP_create_G_from_qdp(g->links);
  QOP_symanzik_1loop_gauge_action(&info, qg, &acts, &actt, coeffs);
  QOP_destroy_G(qg);
  return beta*(acts + actt);
}
#endif

static int
qopqdp_gauge_force(lua_State *L)
{
  qassert(lua_gettop(L)==3||lua_gettop(L)==4);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  force_t *f = qopqdp_force_check(L, 2);
  for(int i=0; i<f->nd; i++) QDP_M_eq_zero(f->force[i], QDP_all);
  QOP_gauge_coeffs_t coeffs;
  get_gauge_coeffs(L, &coeffs, 3);
  double beta = luaL_optnumber(L, 4, 1);
  QOP_info_t info;
  QOP_GaugeField *qg = QOP_create_G_from_qdp(g->links);
  QOP_Force *qf = QOP_create_F_from_qdp(f->force);
  QOP_symanzik_1loop_gauge_force(&info, qg, qf, &coeffs, beta);
  QOP_extract_F_to_qdp(f->force, qf);
  QOP_destroy_G(qg);
  QOP_destroy_F(qf);
#if 0
  QOP_symanzik_1loop_gauge_deriv(&info, qg, qf, &coeffs, beta);
  QOP_extract_F_to_qdp(f->force, qf);
  QOP_destroy_G(qg);
  QOP_destroy_F(qf);
  get_lie_force(f, g);
  get_measure_force(f, g);
#endif
  //check_force(f, g, action, (void*[]){(void*)&coeffs,(void*)&beta});
  f->time = info.final_sec;
  f->flops = info.final_flop;
  return 0;
}

static int
qopqdp_gauge_update(lua_State *L)
{
  qassert(lua_gettop(L)==3);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  force_t *f = qopqdp_force_check(L, 2);
  QLA_Real eps = luaL_checknumber(L, 3);
#if 1
  //get_gauge_links(g);
  QDP_ColorMatrix *m1 = QDP_create_M();
  QDP_ColorMatrix *m2 = QDP_create_M();
  for(int i=0; i<g->nd; i++) {
    QDP_M_eq_r_times_M(m1, &eps, f->force[i], QDP_all);
    QDP_M_eq_exp_M(m2, m1, QDP_all);
    QDP_M_eq_M_times_M(m1, m2, g->links[i], QDP_all);
    QDP_M_eq_M(g->links[i], m1, QDP_all);
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
  if(g->lie) get_gauge_lie(g);
#else
  QLA_Complex ieps;
  QLA_c_eq_r_plus_ir(ieps, 0, -eps);
  for(int i=0; i<g->nd; i++) {
    QDP_M_peq_c_times_M(g->lie[i], &ieps, f->force[i], QDP_all);
    QDP_M_eq_funcia(g->lie[i], make_herm, NULL, QDP_all);
  }
#if 0
  {
    printf0("initial mat:\n");
    QLA_ColorMatrix *cm = QDP_expose_M(g->lie[0]);
    for(int i=0; i<QLA_Nc; i++) {
      for(int j=0; j<QLA_Nc; j++) {
	QLA_Complex *z = &QLA_elem_M(*cm, i, j);
	printf0(" lie%i%i\t%g\t%g\n", i, j, QLA_real(*z), QLA_imag(*z));
      }
    }
    QDP_reset_M(g->lie[0]);
  }
  for(int i=0; i<g->nd; i++) {
    QDP_M_eq_funcia(g->lie[i], make_herm, NULL, QDP_all);
  }
  get_gauge_links(g);
  get_gauge_lie(g);
  {
    printf0("normalized mat:\n");
    QLA_ColorMatrix *cm = QDP_expose_M(g->lie[0]);
    for(int i=0; i<QLA_Nc; i++) {
      for(int j=0; j<QLA_Nc; j++) {
	QLA_Complex *z = &QLA_elem_M(*cm, i, j);
	printf0(" lie%i%i\t%g\t%g\n", i, j, QLA_real(*z), QLA_imag(*z));
      }
    }
    QDP_reset_M(g->lie[0]);
  }
  for(int i=0; i<g->nd; i++) {
    QDP_M_eq_funcia(g->lie[i], make_herm, NULL, QDP_all);
  }
#endif
  get_gauge_links(g);
#endif
  return 0;
}

#if 0
typedef struct {
  QLA_Real eps;
  QLA_ColorMatrix *f;
} UGargs;

static void
gauge_update(QLA_ColorMatrix *g, int i, void *args)
{
  UGargs *a = (UGargs *) args;
  QLA_ColorMatrix m1, m2;
  QLA_M_eq_r_times_M(&m1, &a->eps, &a->f[i]);
  QLA_M_eq_exp_M(&m2, &m1);
  QLA_M_eq_M_times_M(&m1, &m2, g);
  QLA_M_eq_M(g, &m1);
}

static int
qopqdp_gauge_update2(lua_State *L)
{
  qassert(lua_gettop(L)==3);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  force_t *f = qopqdp_force_check(L, 2);
  UGargs a;
  a.eps = luaL_checknumber(L, 3);
  for(int i=0; i<g->nd; i++) {
    a.f = QDP_expose_M(f->force[i]);
    QDP_M_eq_funcia(g->links[i], gauge_update, &a, QDP_all);
    QDP_reset_M(f->force[i]);
  }
  return 0;
}
#endif

// calculate generic loop
static int
qopqdp_gauge_loop(lua_State *L)
{
  qassert(lua_gettop(L)==2);
  gauge_t *g = qopqdp_gauge_check(L, 1);
  int ns; get_table_len(L, 2, &ns);
  int dirs[ns]; get_int_array(L, 2, ns, dirs);

  QLA_Complex z;
  QDP_ColorMatrix *m[2], *temp[2];
  m[0] = QDP_create_M();
  m[1] = QDP_create_M();
  temp[0] = QDP_create_M();
  temp[1] = QDP_create_M();

  int k = 0;
  int intemp = 0;
  QLA_c_eq_r(z, 1);
  QDP_M_eq_c(m[k], &z, QDP_all);
  for(int i=0; i<ns; i++) {
    int d = dirs[i];
    int mu = abs(d) - 1;
    //printf0("i: %i  d: %i  mu: %i\n", i, d, mu);
    if(mu<0 || mu>=g->nd) {
      printf0("%s invalid index[%i] %i (nd=%i)\n", __func__, i, d, g->nd);
      lua_pushnil(L);
      return 1;
    }
    if(d>0) { // shift from backward
      if(intemp) {
	QDP_M_eq_Ma_times_M(m[1-k], g->links[mu], temp[k], QDP_all);
	//QDP_discard_M(temp[k]);
      } else {
	QDP_M_eq_Ma_times_M(m[1-k], g->links[mu], m[k], QDP_all);
      }
      k = 1-k;
      QDP_M_eq_sM(temp[k], m[k], QDP_neighbor[mu], QDP_backward, QDP_all);
      intemp = 1;
    } else { // shift from forward
      if(intemp) {
	QDP_M_eq_M(m[1-k], temp[k], QDP_all);
	//QDP_discard_M(temp[k]);
	k = 1-k;
      }
      //printf("k: %i\n", k);
      //TRACE;
      QDP_M_eq_sM(temp[k], m[k], QDP_neighbor[mu], QDP_forward, QDP_all);
      //TRACE;
      QDP_M_eq_M_times_M(m[1-k], g->links[mu], temp[k], QDP_all);
      //TRACE;
      k = 1-k;
      intemp = 0;
    }
  }
  if(intemp) {
    QDP_M_eq_M(m[1-k], temp[k], QDP_all);
    QDP_discard_M(temp[k]);
    k = 1-k;
  }

  QLA_ColorMatrix cm;
  QDP_m_eq_sum_M(&cm, m[k], QDP_all);
  QLA_C_eq_trace_M(&z, &cm);
  double f = 1.0/(QLA_Nc*QDP_volume());
  lua_pushnumber(L, f*QLA_real(z));
  lua_pushnumber(L, f*QLA_imag(z));

  QDP_destroy_M(m[0]);
  QDP_destroy_M(m[1]);
  QDP_destroy_M(temp[0]);
  QDP_destroy_M(temp[1]);

  return 2;
}

static struct luaL_Reg gauge_reg[] = {
  { "__gc",    qopqdp_gauge_gc },
  { "unit",    qopqdp_gauge_unit },
  { "random",  qopqdp_gauge_random },
  { "load",    qopqdp_gauge_load },
  { "save",    qopqdp_gauge_save },
  { "set",     qopqdp_gauge_set },
  { "copy",    qopqdp_gauge_copy },
  { "checkU",  qopqdp_gauge_checkU },
  { "checkSU", qopqdp_gauge_checkSU },
  { "makeSU",  qopqdp_gauge_makeSU },
  { "action",  qopqdp_gauge_action },
  { "force",   qopqdp_gauge_force },
  { "update",  qopqdp_gauge_update },
  { "loop",    qopqdp_gauge_loop },
  { NULL, NULL}
};

gauge_t *
qopqdp_gauge_create(lua_State* L)
{
  int nd = QDP_ndim();
  gauge_t *g = lua_newuserdata(L, sizeof(gauge_t));
  g->nd = nd;
  g->links = malloc(nd*sizeof(QDP_ColorMatrix*));
  //g->lie = malloc(nd*sizeof(QDP_ColorMatrix*));
  g->lie = NULL;
  for(int i=0; i<nd; i++) {
    g->links[i] = QDP_create_M();
    if(g->lie) g->lie[i] = QDP_create_M();
  }
  if(luaL_newmetatable(L, gmtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, gauge_reg);
  }
  lua_setmetatable(L, -2);
  return g;
}
