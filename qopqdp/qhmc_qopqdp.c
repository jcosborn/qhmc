#include "qhmc_qopqdp_common.h"
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <qmp.h>

QLA_RandomState qopqdp_nrs;
QDP_RandomState *qopqdp_srs = NULL;

void
init_qopqdp(int *argc, char ***argv)
{
  QDP_initialize(argc, argv);
}

void
fini_qopqdp(void)
{
  QDP_finalize();
}

static int
slice_func(int x[], void *args)
{
  int dir = *(int *)args;
  return x[dir];
}

QDP_Subset *
qhmcqdp_get_timeslices(lattice_t *lat)
{
  if(lat->timeslices==NULL) {
    int dir = lat->nd - 1;
    int n = QDP_coord_size(dir);
    lat->timeslices = QDP_create_subset(slice_func,(void*)&dir,sizeof(dir),n);
  }
  return lat->timeslices;
}

QOP_evenodd_t
qopqdp_check_evenodd(lua_State *L, int idx)
{
  lua_pushvalue(L, idx);
  const char *s = luaL_checkstring(L, -1);
  lua_pop(L, 1);
  QOP_evenodd_t eo=QOP_EVENODD;
  switch(s[0]) {
  case 'a': break;
  case 'e': eo = QOP_EVEN; break;
  case 'o': eo = QOP_ODD; break;
  default: qerror(1, "unknown parity %s\n", s);
  }
  return eo;
}

// FIXME: make work for non-default subset
QDP_Subset
qopqdp_opt_subset(lua_State *L, int *idx, int reqd, QDP_Subset def)
{
  lattice_t *lat = qopqdp_get_default_lattice(L);
  lua_pushvalue(L, *idx);
  const char *s;
  if(reqd) {
    s = luaL_checkstring(L, -1);
  } else {
    s = luaL_optstring(L, -1, NULL);
  }
  lua_pop(L, 1);
  QDP_Subset sub = def;
  if(s) {
    switch(s[0]) {
    case 'a': sub = QDP_all_L(lat->qlat); break;
    case 'e': sub = QDP_even_L(lat->qlat); break;
    case 'o': sub = QDP_odd_L(lat->qlat); break;
    case 't': 
      if(strncmp(s,"timeslice",9)==0) {
	int t;
	int n = sscanf(s+9,"%i",&t);
	if(n && t>=0 && t<QDP_coord_size_L(lat->qlat,QDP_ndim_L(lat->qlat)-1))
	  sub = qhmcqdp_get_timeslices(lat)[t];
      }
      break;
    }
    if(sub==NULL) qlerror0(L, 1, "unknown subset %s\n", s);
    (*idx)++;
  }
  return sub;
}

static int
qopqdp_master(lua_State* L)
{
  qassert(lua_gettop(L)==0);
  lua_pushboolean(L, QDP_this_node==0);
  return 1;
}

static int
qopqdp_dtime(lua_State* L)
{
  qassert(lua_gettop(L)==0);
  lua_pushnumber(L, QMP_time());
  return 1;
}

static void
qopqdp_set_default_lattice(lua_State *L, int idx)
{
  lua_pushvalue(L, idx);
  lua_setfield(L, LUA_REGISTRYINDEX, "defaultLattice");
}

lattice_t *
qopqdp_get_default_lattice(lua_State *L)
{
  lua_getfield(L, LUA_REGISTRYINDEX, "defaultLattice");
  lattice_t *lat = qopqdp_opt_lattice(L, (int[]){-1}, 0, NULL);
  lua_pop(L, 1);
  return lat;
}

static int
qopqdp_lattice(lua_State *L)
{
  int nargs = lua_gettop(L);
  if(nargs==0) {
    int nd = QDP_ndim();
    int lat[nd];
    QDP_latsize(lat);
    push_int_array(L, nd, lat);
  } else {
    int nd;
    get_table_len(L, -1, &nd);
    int size[nd];
    get_int_array(L, -1, nd, size);
    lattice_t *lat = qopqdp_lattice_create(L, nd, size);
    lua_pushvalue(L, -1);
    lat->ref = luaL_ref(L, LUA_REGISTRYINDEX); // prevent gc
    if(QDP_get_default_lattice()==NULL) { // no default lattice
      qopqdp_set_default_lattice(L, -1);
      QDP_set_default_lattice(lat->qlat);
      qopqdp_srs = QDP_create_S();
      QLA_use_milc_gaussian = 1;
      QOP_layout_t qoplayout;
      qoplayout.latdim = nd;
      qoplayout.latsize = size;
      qoplayout.machdim = -1;
      QOP_init(&qoplayout);
    }
  }
  return 1;
}

static int
qopqdp_profile(lua_State* L)
{
  int nargs = lua_gettop(L);
  qassert(nargs>=0 || nargs<=1);
  int r;
  if(nargs==0) {
    r = QDP_profcontrol(0);
    QDP_profcontrol(r);
  } else {
    int v = luaL_checkinteger(L, 1);
    r = QDP_profcontrol(v);
  }
  lua_pushinteger(L, r);
  return 1;
}

static int
qopqdp_verbosity(lua_State* L)
{
  int nargs = lua_gettop(L);
  qassert(nargs>=0 || nargs<=1);
  int r;
  if(nargs==0) {
    r = QOP_verbose(0);
    QOP_verbose(r);
  } else {
    int v = luaL_checkinteger(L, 1);
    r = QOP_verbose(v);
  }
  lua_pushinteger(L, r);
  return 1;
}

static int
qopqdp_blocksize(lua_State* L)
{
  BEGIN_ARGS;
  OPT_INT(new,0);
  END_ARGS;
  int old = QDP_get_block_size();
  lua_pushinteger(L, old);
  if(new>0) QDP_set_block_size(new);
  return 1;
}

static int
qopqdp_readGroupSize(lua_State* L)
{
  BEGIN_ARGS;
  OPT_INT(new,0);
  END_ARGS;
  int old = QDP_set_read_group_size(new);
  lua_pushinteger(L, old);
  if(new<=0) QDP_set_read_group_size(old);
  return 1;
}

static int
qopqdp_writeGroupSize(lua_State* L)
{
  BEGIN_ARGS;
  OPT_INT(new,0);
  END_ARGS;
  int old = QDP_set_write_group_size(new);
  lua_pushinteger(L, old);
  if(new<=0) QDP_set_write_group_size(old);
  return 1;
}

static void
lex_int(QLA_Int *li, int coords[])
{
  int nd = QDP_ndim();
  int t = coords[nd-1];
  for(int i=nd-2; i>=0; i--) {
    t = t*QDP_coord_size(i) + coords[i];
  }
  *li = t;
}

static int
qopqdp_seed(lua_State* L)
{
  qassert(lua_gettop(L)==1);
  int seed = luaL_checkint(L, 1);
  QDP_Int *li = QDP_create_I();
#ifdef QHMC_REPRO_UNIFORM
  QDP_I_eq_i(li, (int[1]){0x55555555}, QDP_all);
#else
  QDP_I_eq_func(li, lex_int, QDP_all);
#endif
  QDP_S_eq_seed_i_I(qopqdp_srs, seed, li, QDP_all);
  QDP_destroy_I(li);
  QLA_Int i = QDP_volume()+QDP_this_node;
  QLA_S_eq_seed_i_I(&qopqdp_nrs, seed, &i);
  return 0;
}

static int
qopqdp_random(lua_State *L)
{
  qassert(lua_gettop(L)==0);
  double r=0;
  if(QDP_this_node==0) {
    QLA_Real t;
    QLA_R_eq_random_S(&t, &qopqdp_nrs);
    r = t;
  }
  //QMP_sum_double(&r);
  QMP_broadcast(&r, sizeof(r));
  lua_pushnumber(L, r);
  return 1;
}

static int
qopqdp_cscalar(lua_State *L)
{
  BEGIN_ARGS;
  OPT_LATTICE(lat, NULL);
  END_ARGS;
  qopqdp_cscalar_create(L, lat);
  return 1;
}

// 1: precision
// 2: nc
// 3: lattice
static int
qopqdp_gauge(lua_State *L)
{
  BEGIN_ARGS;
  OPT_STRING(precision, "D");
  OPT_INT(nc, 0);
  OPT_LATTICE(lat, NULL);
  END_ARGS;
  if(*precision=='F') {
    qopqdp_gaugeF_create(L, nc, lat);
  } else {
    qopqdp_gaugeD_create(L, nc, lat);
  }
  return 1;
}

// 1: precision
// 2: nc
// 3: lattice
static int
qopqdp_force(lua_State *L)
{
  BEGIN_ARGS;
  OPT_STRING(precision, "D");
  OPT_INT(nc, 0);
  OPT_LATTICE(lat, NULL);
  END_ARGS;
  if(*precision=='F') {
    qopqdp_forceF_create(L, nc, lat);
  } else {
    qopqdp_forceD_create(L, nc, lat);
  }
  return 1;
}

static int
qopqdp_asqtad(lua_State *L)
{
  qassert(lua_gettop(L)==0);
  qopqdp_asqtad_create(L);
  return 1;
}

static int
qopqdp_hisq(lua_State *L)
{
  qassert(lua_gettop(L)==0);
  qopqdp_hisq_create(L);
  return 1;
}

static int
qopqdp_wilson(lua_State *L)
{
  qassert(lua_gettop(L)==0);
  qopqdp_wilson_create(L);
  return 1;
}

static int
qopqdp_dw(lua_State *L)
{
  qassert(lua_gettop(L)==0);
  qopqdp_dw_create(L);
  return 1;
}

// 1: filename
// 2: (opt) lattice
// return: reader, file metadata
static int
qopqdp_reader(lua_State *L)
{
  BEGIN_ARGS;
  GET_STRING(fn);
  OPT_LATTICE(lat, NULL);
  END_ARGS;
  qopqdp_reader_create(L, fn, lat);
  return 2;
}

// 1: filename
// 2: metadata
static int
qopqdp_writer(lua_State *L)
{
  int narg = lua_gettop(L);
  qassert(narg==2);
  lua_pushvalue(L, 1);
  const char *fn = luaL_checkstring(L, -1);
  lua_pop(L, 1);
  lua_pushvalue(L, 2);
  const char *md = luaL_checkstring(L, -1);
  lua_pop(L, 1);
  qopqdp_writer_create(L, fn, md);
  return 1;
}

static int
qopqdp_remapout(lua_State *L)
{
  qassert(lua_gettop(L)==1);
  lua_pushvalue(L, -1);
  const char *s = luaL_checkstring(L, -1);
  lua_pop(L, 1);
  int fd = creat(s, 0666);
  fflush(stdout);
  dup2(fd, 1);
  fflush(stderr);
  dup2(fd, 2);
  close(fd);
  return 0;
}

static struct luaL_Reg qopqdp_reg[] = {
  { "master",         qopqdp_master },
  { "dtime",          qopqdp_dtime },
  { "lattice",        qopqdp_lattice },
  { "profile",        qopqdp_profile },
  { "verbosity",      qopqdp_verbosity },
  { "blocksize",      qopqdp_blocksize },
  { "readGroupSize",  qopqdp_readGroupSize },
  { "writeGroupSize", qopqdp_writeGroupSize },
  { "seed",           qopqdp_seed },
  { "random",         qopqdp_random },
  { "cscalar",        qopqdp_cscalar },
  { "gauge",          qopqdp_gauge },
  { "force",          qopqdp_force },
  { "asqtad",         qopqdp_asqtad },
  { "hisq",           qopqdp_hisq },
  { "wilson",         qopqdp_wilson },
  { "dw",             qopqdp_dw },
  { "reader",         qopqdp_reader },
  { "writer",         qopqdp_writer },
  { "remapout",       qopqdp_remapout },
  { NULL, NULL}
};

void
open_qopqdp(lua_State* L)
{
  luaL_register(L, "qopqdp", qopqdp_reg);
  int jobnum = QMP_get_job_number();
  int numjobs = QMP_get_number_of_jobs();
  lua_pushinteger(L, jobnum);
  lua_setglobal(L, "jobnum");
  lua_pushinteger(L, numjobs);
  lua_setglobal(L, "numjobs");
  lua_getglobal(L, "qopqdp");
#ifdef DEFAULTNC
  lua_pushinteger(L, DEFAULTNC);
  lua_setfield(L, -2, "Nc");
#endif
  lua_pop(L, 1);
  QDP_set_read_group_size(64);
  QDP_set_write_group_size(64);
  open_qopqdp_smear(L);
}
