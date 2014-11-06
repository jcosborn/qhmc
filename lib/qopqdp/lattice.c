#include "qhmc_qopqdp_common.h"

static char *mtname = "qopqdp.lattice";

lattice_t *
qopqdp_opt_lattice(lua_State *L, int *idx, int required, lattice_t *def)
{
  lattice_t *lat;
  if(required) {
    lat = luaL_checkudata(L, *idx, mtname);
    (*idx)++;
  } else {
    lat = luaL_testudata(L, *idx, mtname);
    if(lat==NULL) lat = def;
    else (*idx)++;
  }
  return lat;
}

static int
qopqdp_lattice_gc(lua_State *L)
{
  lattice_t *lat = qopqdp_check_lattice(L, -1);
  if(lat->timeslices) {
    QDP_destroy_subset(lat->timeslices);
    lat->timeslices = NULL;
  }
  if(lat->staggered) {
    QDP_destroy_subset(lat->staggered);
    lat->staggered = NULL;
  }
  for(int i=0; i<=lat->nd; i++) {
    if(lat->eodir[i]) {
      QDP_destroy_subset(lat->eodir[i]);
      lat->eodir[i] = NULL;
    }
  }
  free(lat->eodir);
  //if(lat->rs) {
  //QDP_destroy_S(lat->rs);
  //lat->rs = NULL;
  //}
  if(lat->qlat!=QDP_get_default_lattice()) {
    QDP_destroy_lattice(lat->qlat);
    lat->qlat = NULL;
  }
  return 0;
}

static int
qopqdp_lattice_len(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(l);
  nextarg++; // ignore second argument
  END_ARGS;
  lua_pushinteger(L, l->nd);
  return 1;
}

static int
qopqdp_lattice_call(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(l);
  OPT_INT(dim, 0);
  END_ARGS;
  if(dim>0) {
    int s = QDP_coord_size_L(l->qlat, dim-1);
    lua_pushinteger(L, s);
  } else {
    int nd = QDP_ndim_L(l->qlat);
    int x[nd];
    QDP_latsize_L(l->qlat, x);
    qhmc_push_int_array(L, nd, x);
  }
  return 1;
}

static int
qopqdp_lattice_volume(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(l);
  END_ARGS;
  int v = QDP_volume_L(l->qlat);
  lua_pushinteger(L, v);
  return 1;
}

static int
qopqdp_lattice_seed(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(l);
  GET_INT(seed);
  OPT_INT(uniform, -1);
  OPT_SUBSET(sub, l, QDP_all_L(l->qlat));
  END_ARGS;
  if(l->rs==NULL) {
    qopqdp_rstate_t *r = qopqdp_rstate_create(L, l);
    l->rs = r->field;
    QHMC_USERTABLE_SETFIELD(L, 1, "rs");
  }
  qhmc_qopqdp_seed_func(l->rs, seed, uniform, sub);
  return 0;
}

static int
qopqdp_lattice_get_rstate(lua_State *L)
{
  //BEGIN_ARGS;
  //GET_LATTICE(l);
  //END_ARGS;
  QHMC_USERTABLE_GETFIELD(L, 1, "rs");
  return 1;
}

static int
qopqdp_lattice_set_rstate(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(l);
  GET_QOPQDP_RSTATE(rs);
  END_ARGS;
  l->rs = rs->field;
  QHMC_USERTABLE_SETFIELD(L, 1, "rs");
  return 0;
}

// 1: lattice
// 2: string or function
// return: subset or array of subsets
static int
qopqdp_lattice_subset(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  OPT_STRING(sn, NULL);
  OPT_FUNCTION_INDEX(fidx, 0);
  END_ARGS;
  int n = 0;
  QDP_Subset *qs = NULL;
  subsetGroup_t *group = NULL;
  if(fidx==0) {
    qs = qhmc_qopqdp_qsubset_from_string(L, lat, sn, &n);
  } else {
    // create qdp subsets
    // create subsetGroup
  }
  if(n==0) {
    lua_pushnil(L);
  } else if(n==1) {
    qopqdp_subset_create(L, qs[0], group);
  } else {
    lua_createtable(L, n, 0);
    for(int i=0; i<n; i++) {
      qopqdp_subset_create(L, qs[i], group);
      lua_rawseti(L, -2, i+1);
    }
  }
  return 1;
}

// 1: lattice
// 2: filename
// return: reader, file metadata
static int
qopqdp_lattice_reader(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  GET_STRING(fn);
  END_ARGS;
  qopqdp_reader_create(L, fn, lat);
  return 2;
}

// 1: lattice
// 2: filename
// 3: metadata
static int
qopqdp_lattice_writer(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  GET_STRING(fn);
  GET_STRING(md);
  END_ARGS;
  qopqdp_writer_create(L, fn, md, lat);
  return 1;
}

static int
qopqdp_rstate(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  END_ARGS;
  qopqdp_rstate_create(L, lat);
  return 1;
}

static int
qopqdp_real(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  OPT_STRING(precision, lat->defaultPrecision);
  END_ARGS;
  if(*precision=='F') {
    qopqdp_realF_create(L, lat);
  } else {
    qopqdp_realD_create(L, lat);
  }
  return 1;
}

static int
qopqdp_complex(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  OPT_STRING(precision, lat->defaultPrecision);
  END_ARGS;
  if(*precision=='F') {
    qopqdp_complexF_create(L, lat);
  } else {
    qopqdp_complexD_create(L, lat);
  }
  return 1;
}

static int
qopqdp_colorVector(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  OPT_STRING(precision, lat->defaultPrecision);
  OPT_INT(nc, lat->defaultNc);
  END_ARGS;
  if(*precision=='F') {
    qopqdp_cvectorF_create(L, nc, lat);
  } else {
    qopqdp_cvectorD_create(L, nc, lat);
  }
  return 1;
}

static int
qopqdp_diracFermion(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  OPT_STRING(precision, lat->defaultPrecision);
  OPT_INT(nc, lat->defaultNc);
  END_ARGS;
  if(*precision=='F') {
    qopqdp_dfermionF_create(L, nc, lat);
  } else {
    qopqdp_dfermionD_create(L, nc, lat);
  }
  return 1;
}

static int
qopqdp_colorMatrix(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  OPT_STRING(precision, lat->defaultPrecision);
  OPT_INT(nc, lat->defaultNc);
  END_ARGS;
  if(*precision=='F') {
    qopqdp_cmatrixF_create(L, nc, lat);
  } else {
    qopqdp_cmatrixD_create(L, nc, lat);
  }
  return 1;
}

static int
qopqdp_gauge(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  OPT_STRING(precision, lat->defaultPrecision);
  OPT_INT(nc, lat->defaultNc);
  END_ARGS;
  if(*precision=='F') {
    qopqdp_gaugeF_create(L, nc, lat);
  } else {
    qopqdp_gaugeD_create(L, nc, lat);
  }
  return 1;
}

static int
qopqdp_force(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  OPT_STRING(precision, lat->defaultPrecision);
  OPT_INT(nc, lat->defaultNc);
  END_ARGS;
  if(*precision=='F') {
    //qopqdp_forceF_create(L, nc, lat);
    qopqdp_gaugeF_create(L, nc, lat);
  } else {
    //qopqdp_forceD_create(L, nc, lat);
    qopqdp_gaugeD_create(L, nc, lat);
  }
  return 1;
}

static int
qopqdp_asqtad(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  OPT_INT(nc, lat->defaultNc);
  END_ARGS;
  qopqdp_asqtad_create(L, nc, lat);
  return 1;
}

static int
qopqdp_hisq(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  OPT_INT(nc, lat->defaultNc);
  END_ARGS;
  qopqdp_hisq_create(L, nc, lat);
  return 1;
}

static int
qopqdp_wilson(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  OPT_INT(nc, lat->defaultNc);
  END_ARGS;
  qopqdp_wilson_create(L, nc, lat);
  return 1;
}

static int
qopqdp_dw(lua_State *L)
{
  BEGIN_ARGS;
  GET_LATTICE(lat);
  OPT_INT(nc, lat->defaultNc);
  END_ARGS;
  qopqdp_dw_create(L, nc, lat);
  return 1;
}

static struct luaL_Reg lattice_reg[] = {
  { "__gc",           qopqdp_lattice_gc },
  { "__len",          qopqdp_lattice_len },
  { "__call",         qopqdp_lattice_call },
  { "volume",         qopqdp_lattice_volume },
  { "seed",           qopqdp_lattice_seed },
  { "subset",         qopqdp_lattice_subset },
  { "reader",         qopqdp_lattice_reader },
  { "writer",         qopqdp_lattice_writer },
  { "getRstate",      qopqdp_lattice_get_rstate },
  { "setRstate",      qopqdp_lattice_set_rstate },
  { "rstate",         qopqdp_rstate },
  { "real",           qopqdp_real },
  { "complex",        qopqdp_complex },
  { "colorVector",    qopqdp_colorVector },
  { "diracFermion",   qopqdp_diracFermion },
  { "colorMatrix",    qopqdp_colorMatrix },
  { "gauge",          qopqdp_gauge },
  { "force",          qopqdp_force },
  { "asqtad",         qopqdp_asqtad },
  { "hisq",           qopqdp_hisq },
  { "wilson",         qopqdp_wilson },
  { "dw",             qopqdp_dw },
  { NULL, NULL}
};

// subsets: 
// int ndirs, dirs[ndirs], origin[ndirs], block[ndirs]
// int ndigits, base[ndigits]
// double digitcoeffs[ndigits]
// * k = dirs[i]  (i=1..ndirs)
// * y[i] = (x[k]-origin[i])//block[i]
// * digit[j] = floor(sum_l dc[j][l]*y[l]) % base[j]

lattice_t *
qopqdp_lattice_wrap(lua_State *L, QDP_Lattice *qlat,
		    char *defPrec, int defNc, int doGC)
{
  int nd = QDP_ndim_L(qlat);
  lattice_t *lat = lua_newuserdata(L, sizeof(lattice_t));
  QHMC_USERTABLE_CREATE(L, -1);
  lat->qlat = qlat;
  lat->timeslices = NULL;
  lat->staggered = NULL;
  lat->eodir = malloc((nd+1)*sizeof(QDP_Subset *));
  for(int i=0; i<=nd; i++) lat->eodir[i] = NULL;
  lat->rs = NULL;
  lat->nd = nd;
  lat->ref = -1;
  lat->defaultPrecision = defPrec;
  lat->defaultNc = defNc;
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, lattice_reg);
  }
  lua_setmetatable(L, -2);
  return lat;
}

lattice_t *
qopqdp_create_lattice(lua_State *L, int nd, int size[],
		      char *defPrec, int defNc)
{
  QDP_Lattice *qlat = QDP_create_lattice(NULL, NULL, nd, size);
  return qopqdp_lattice_wrap(L, qlat, defPrec, defNc, 1);
}
