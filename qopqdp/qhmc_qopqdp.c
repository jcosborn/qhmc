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

void
get_bool_array(lua_State *L, int idx, int n, int *a)
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qassert(n==lua_objlen(L, idx));
  for(int i=0; i<n; i++) {
    lua_rawgeti(L, idx, i+1);
    a[i] = lua_toboolean(L, -1);
    lua_pop(L, 1);
  }
}

void
get_int_array(lua_State *L, int idx, int n, int *a)
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qassert(n==lua_objlen(L, idx));
  for(int i=0; i<n; i++) {
    lua_rawgeti(L, idx, i+1);
    a[i] = lua_tointeger(L, -1);
    lua_pop(L, 1);
  }
}

void
push_int_array(lua_State *L, int n, int *a)
{
  lua_createtable(L, n, 0);
  for(int i=0; i<n; i++) {
    lua_pushinteger(L, a[i]);
    lua_rawseti(L, -2, i+1);
  }
}

void
get_float_array(lua_State *L, int idx, int n, float *a)
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qassert(n==lua_objlen(L, idx));
  for(int i=0; i<n; i++) {
    lua_rawgeti(L, idx, i+1);
    a[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
}

void
get_double_array(lua_State *L, int idx, int n, double *a)
{
  luaL_checktype(L, idx, LUA_TTABLE);
  qassert(n==lua_objlen(L, idx));
  for(int i=0; i<n; i++) {
    lua_rawgeti(L, idx, i+1);
    a[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
}

void
push_double_array(lua_State *L, int n, double *a)
{
  lua_createtable(L, n, 0);
  for(int i=0; i<n; i++) {
    lua_pushnumber(L, a[i]);
    lua_rawseti(L, -2, i+1);
  }
}

const char *
tableGetString(lua_State *L, int idx, char *key)
{
  lua_getfield(L, idx, key);
  const char *s = luaL_checkstring(L, -1);
  lua_pop(L, 1);
  return s;
}

int
slice_func(int x[], void *args)
{
  int dir = *(int *)args;
  return x[dir];
}

QDP_Subset *qhmcqdp_timeslices = NULL;
QDP_Subset *

qhmcqdp_get_timeslices(void)
{
  if(qhmcqdp_timeslices==NULL) {
    int dir = 3;
    int n = QDP_coord_size(dir);
    qhmcqdp_timeslices = QDP_create_subset(slice_func, (void*)&dir, sizeof(dir), n);
  }
  return qhmcqdp_timeslices;
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
  default: qerror("unknown parity %s\n", s);
  }
  return eo;
}

QDP_Subset
qopqdp_check_subset(lua_State *L, int idx)
{
  lua_pushvalue(L, idx);
  const char *s = luaL_checkstring(L, -1);
  lua_pop(L, 1);
  QDP_Subset sub = QDP_all;
  switch(s[0]) {
  case 'a': break;
  case 'e': sub = QDP_even; break;
  case 'o': sub = QDP_odd; break;
  default: qerror("unknown subset %s\n", s);
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

static int
qopqdp_lattice(lua_State *L)
{
  int nret = 0;
  int nargs = lua_gettop(L);
  if(nargs==0) {
    int nd = QDP_ndim();
    int lat[nd];
    QDP_latsize(lat);
    push_int_array(L, nd, lat);
    nret = 1;
  } else {
    int nd;
    get_table_len(L, -1, &nd);
    int lat[nd];
    get_int_array(L, -1, nd, lat);
    QDP_set_latsize(nd, lat);
    QDP_create_layout();

    qopqdp_srs = QDP_create_S();
    QLA_use_milc_gaussian = 1;

    QOP_layout_t qoplayout;
    qoplayout.latdim = nd;
    qoplayout.latsize = lat;
    qoplayout.machdim = -1;
    QOP_init(&qoplayout);
  }
  return nret;
}

static int
qopqdp_profile(lua_State* L)
{
  qassert(lua_gettop(L)==1);
  int v = luaL_checkinteger(L, 1);
  QDP_profcontrol(v);
  return 0;
}

static int
qopqdp_verbosity(lua_State* L)
{
  qassert(lua_gettop(L)==1);
  int v = luaL_checkinteger(L, 1);
  QOP_verbose(v);
  return 0;
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
qopqdp_random(lua_State* L)
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
qopqdp_gauge(lua_State* L)
{
  qassert(lua_gettop(L)==0);
  qopqdp_gauge_create(L);
  return 1;
}

static int
qopqdp_force(lua_State* L)
{
  qassert(lua_gettop(L)==0);
  qopqdp_force_create(L);
  return 1;
}

static int
qopqdp_asqtad(lua_State* L)
{
  qassert(lua_gettop(L)==0);
  qopqdp_asqtad_create(L);
  return 1;
}

static int
qopqdp_hisq(lua_State* L)
{
  qassert(lua_gettop(L)==0);
  qopqdp_hisq_create(L);
  return 1;
}

static int
qopqdp_wilson(lua_State* L)
{
  qassert(lua_gettop(L)==0);
  qopqdp_wilson_create(L);
  return 1;
}

static int
qopqdp_dw(lua_State* L)
{
  qassert(lua_gettop(L)==0);
  qopqdp_dw_create(L);
  return 1;
}

static int
qopqdp_remapout(lua_State* L)
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
  { "master",    qopqdp_master },
  { "dtime",     qopqdp_dtime },
  { "lattice",   qopqdp_lattice },
  { "profile",   qopqdp_profile },
  { "verbosity", qopqdp_verbosity },
  { "seed",      qopqdp_seed },
  { "random",    qopqdp_random },
  { "gauge",     qopqdp_gauge },
  { "force",     qopqdp_force },
  { "asqtad",    qopqdp_asqtad },
  { "hisq",      qopqdp_hisq },
  { "wilson",    qopqdp_wilson },
  { "dw",        qopqdp_dw },
  { "remapout",  qopqdp_remapout },
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
  lua_pushinteger(L, QLA_Nc);
  lua_setfield(L, -2, "Nc");
  lua_pop(L, 1);
  open_qopqdp_smear(L);
}
