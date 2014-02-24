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
  if(lat->rs) {
    QDP_destroy_S(lat->rs);
    lat->rs = NULL;
  }
  if(lat->qlat!=QDP_get_default_lattice()) {
    QDP_destroy_lattice(lat->qlat);
    lat->qlat = NULL;
  }
  return 0;
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
    l->rs = QDP_create_S_L(l->qlat);
  }
  qhmc_qopqdp_seed_func(l->rs, seed, uniform, sub);
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
qopqdp_cmatrix(lua_State *L)
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

static struct luaL_Reg lattice_reg[] = {
  { "__gc",           qopqdp_lattice_gc },
  { "seed",           qopqdp_lattice_seed },
  { "subset",         qopqdp_lattice_subset },
  { "reader",         qopqdp_lattice_reader },
  { "writer",         qopqdp_lattice_writer },
  { "rstate",         qopqdp_rstate },
  { "real",           qopqdp_real },
  { "complex",        qopqdp_complex },
  { "cmatrix",        qopqdp_cmatrix },
  { NULL, NULL}
};

lattice_t *
qopqdp_create_lattice(lua_State *L, int nd, int size[])
{
  lattice_t *lat = lua_newuserdata(L, sizeof(lattice_t));
  lat->qlat = QDP_create_lattice(NULL, NULL, nd, size);
  lat->timeslices = NULL;
  lat->staggered = NULL;
  lat->eodir = malloc((nd+1)*sizeof(QDP_Subset *));
  for(int i=0; i<=nd; i++) lat->eodir[i] = NULL;
  lat->rs = NULL;
  lat->nd = nd;
  lat->ref = -1;
  lat->defaultPrecision = "D";
  lat->defaultNc = QOPQDP_DEFAULTNC;
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, lattice_reg);
  }
  lua_setmetatable(L, -2);
  return lat;
}
