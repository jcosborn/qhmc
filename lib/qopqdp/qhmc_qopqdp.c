/// Base routines for QOPQDP library support.
// @module qopqdp
#include "qhmc_qopqdp_common.h"
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <qmp.h>

QLA_RandomState qopqdp_nrs;
QDP_RandomState *qopqdp_srs = NULL;

#ifdef QOPQDP_DEFAULTNC
static int defaultNc = QOPQDP_DEFAULTNC;
#else
static int defaultNc = 3;
#endif

void
qhmc_init_qopqdp(int *argc, char ***argv)
{
  QDP_initialize(argc, argv);
  QDP_profcontrol(0);
}

void
qhmc_fini_qopqdp(void)
{
  QDP_finalize();
}

QOP_evenodd_t
qopqdp_opt_evenodd(lua_State *L, int *idx, int required, QOP_evenodd_t def)
{
  const char *s = "(NULL)";
  int valid = 0;
  QOP_evenodd_t eo = def;
  lua_pushvalue(L, *idx);
  if(lua_type(L,-1)==LUA_TSTRING) {
    s = luaL_checkstring(L, -1);
    valid = 1;
    switch(s[0]) {
    case 'a': break;
    case 'e': eo = QOP_EVEN; break;
    case 'o': eo = QOP_ODD; break;
    default: valid = 0;
    }
    if(valid) (*idx)++;
  }
  lua_pop(L, 1);
  if(required && !valid) {
    qlerror(L, 1, "unknown parity %s\n", s);
  }
  return eo;
}

static int
subCopyHyper(QDP_Lattice *rlat, int x[], void *args)
{
  int color = 0;
  int nd = QDP_ndim_L(rlat);
  int *sublen = (int *)args;
  int *rof = sublen + nd;
  int *rls = rof + nd;
  for(int i=0; i<nd; i++) {
    int k = (x[i] - rof[i] + rls[i])%rls[i];
    if(k>=sublen[i]) color = 1;
  }
#if 0
  if(color==0) {
    printf("SUB:");
    for(int i=0; i<nd; i++) printf(" %i", x[i]);
    printf("\n");
  }
#endif
  return color;
}

// shift map from slat into rlat with offsets passed in args
static void
mapCopyHyper(QDP_Lattice *rlat, QDP_Lattice *slat, int rx[], int sx[],
             int *num, int idx, QDP_ShiftDir fb, void *args)
{
  int rnd = QDP_ndim_L(rlat);
  int snd = QDP_ndim_L(slat);
  *num = 1;
  if(fb==QDP_forward) {
    int *rof = (int *)args;
    int *rls = rof + rnd;
    int *sd = rls + rnd;
    int *sof = sd + rnd;
    int *sls = sof + snd;
#if 0
    {
      printf("rof:");
      for(int i=0; i<rnd; i++) printf(" %i", rof[i]);
      printf("rls:");
      for(int i=0; i<rnd; i++) printf(" %i", rls[i]);
      printf("sd: ");
      for(int i=0; i<rnd; i++) printf(" %i", sd[i]);
      printf("sof:");
      for(int i=0; i<snd; i++) printf(" %i", sof[i]);
      printf("sls:");
      for(int i=0; i<snd; i++) printf(" %i", sls[i]);
    }
#endif
    for(int j=0; j<snd; j++) sx[j] = sof[j];
    for(int i=0; i<rnd; i++) {
      int k = (rx[i] - rof[i] + rls[i])%rls[i];
      int j = sd[i];
      sx[j] = (k + sx[j])%sls[j];
      if(k>=sls[j]) *num = 0;
    }
#if 0
    if(*num) {
      printf("FWD:");
      for(int i=0; i<rnd; i++) printf(" %i", rx[i]);
      printf(" <-");
      for(int i=0; i<snd; i++) printf(" %i", sx[i]);
      printf("\n");
    }
#endif
  } else { // QDP_backward
    int *sof = (int *)args;
    int *sls = sof + snd;
    int *sd = sls + snd;
    int *rof = sd + snd;
    int *rls = rof + rnd;
    for(int j=0; j<snd; j++) {
      int i = sd[j];
      int k = 0;
      if(i>=0&&i<rnd) k = (rx[i] - rof[i] + rls[i])%rls[i];
      sx[j] = (k + sof[j])%sls[j];
      if(k>=sls[j]) *num = 0;
    }
#if 0
    if(*num) {
      printf("BCK:");
      for(int i=0; i<snd; i++) printf(" %i", sx[i]);
      printf(" ->");
      for(int i=0; i<rnd; i++) printf(" %i", rx[i]);
      printf("\n");
    }
#endif
  }
}

// creates shifts and maps for copying hypercubic region of size rlen
// from slat to rlat.
// the point soff in slat will get copied to roff in rlat.
// subsequent points in rlat will correspond to the permuted directions
// in slat given by sdir
// j = sdir[i]
// r[i] = ( roff[i] + (s[j]-soff[j]+ss[j])%ss[j] )%rs[i]
// s[j] = ( soff[j] + (r[i]-roff[i]+rs[i])%rs[i] )%ss[j]
void
qhmc_qopqdp_getCopyHyper(QDP_Shift *map, QDP_Subset **subset,
			 QDP_Lattice *rlat, int roff[], int rlen[], int sdir[],
			 QDP_Lattice *slat, int soff[], int num)
{
  int rnd = QDP_ndim_L(rlat);
  int snd = QDP_ndim_L(slat);
  //int sublen[rnd], rof[rnd], rs[rnd], sd[rnd], sof[snd], ss[snd];
  int sublen[4*rnd+2*snd], *rof, *rs, *sd, *sof, *ss, nsub[rnd], nsubs=1;
  rof = sublen + rnd;
  rs = rof + rnd;
  sd = rs + rnd;
  sof = sd + rnd;
  ss = sof + snd;
  QDP_latsize_L(rlat, rs);
  QDP_latsize_L(slat, ss);
  // get subvolume size
  for(int i=0; i<rnd; i++) {
    sd[i] = sdir[i];
    int j = sdir[i];
    sublen[i] = rlen[i];
    if(sublen[i]>rs[i]) sublen[i] = rs[i];
    if(j<0||j>=snd) sublen[i] = 1;
    else if(sublen[i]>ss[j]) sublen[i] = ss[j];
    nsub[i] = (rlen[i]+sublen[i]-1)/sublen[i];
    nsubs *= nsub[i];
  }
  if(num<0||num>=nsubs) {
    *map = NULL;
    *subset = NULL;
    return;
  }
  // calculate which subvolume we will work on and adjust size if necessary
  int n = num;
  for(int j=0; j<snd; j++) sof[j] = soff[j];
  for(int i=0; i<rnd; i++) {
    int j = sdir[i];
    int k = n%nsub[i];
    n = n/nsub[i];
    int off = k*sublen[i];
    rof[i] = (roff[i]+off)%rs[i];
    if(j>=0&&j<snd) sof[j] = (sof[j]+off)%ss[j];
    if(sublen[i]>(rlen[i]-off)) sublen[i] = rlen[i]-off;
  }
  *subset=QDP_create_subset_L(rlat,subCopyHyper,sublen,3*rnd*sizeof(int),2);
  *map=QDP_create_map_L(rlat,slat,mapCopyHyper,rof,(3*rnd+2*snd)*sizeof(int));
}

int
qhmc_qopqdp_master(void)
{
  return QDP_this_node==0;
}

static int
qopqdp_master(lua_State *L)
{
  BEGIN_ARGS;
  END_ARGS;
  lua_pushboolean(L, qhmc_qopqdp_master());
  return 1;
}

static int
qopqdp_rank(lua_State *L)
{
  BEGIN_ARGS;
  END_ARGS;
  lua_pushinteger(L, QDP_this_node);
  return 1;
}

static int
qopqdp_dtime(lua_State *L)
{
  BEGIN_ARGS;
  END_ARGS;
  lua_pushnumber(L, QMP_time());
  return 1;
}

static int
qopqdp_defaultNc(lua_State *L)
{
  BEGIN_ARGS;
  OPT_INT(nc, 0);
  END_ARGS;
  if(nc>0) {
#if QOP_Colors != 'N'
    qlassert(L, nc==QOP_Colors);
#endif
    defaultNc = nc;
    return 0;
  }
  lua_pushinteger(L, defaultNc);
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
qopqdp_defaultLattice(lua_State *L)
{
  BEGIN_ARGS;
  OPT_LATTICE(l, NULL);
  END_ARGS;
  if(l) {
    qopqdp_set_default_lattice(L, 1);
    return 0;
  }
  lua_getfield(L, LUA_REGISTRYINDEX, "defaultLattice");
  return 1;
}

/// Create lattice object (or return size of default lattice).
//  @function lattice
//  @param[opt] dimensions array of lattice dimensions
//  @return lattice object if dimensions given, dimensions of default lattice
//  otherwise
static int
qopqdp_lattice(lua_State *L)
{
  int nargs = lua_gettop(L);
  if(nargs==0) {
    int nd = QDP_ndim();
    int lat[nd];
    QDP_latsize(lat);
    qhmc_push_int_array(L, nd, lat);
  } else {
    int nd;
    get_table_len(L, -1, &nd);
    int size[nd];
    qhmc_get_int_array(L, -1, nd, size);
    lattice_t *lat = qopqdp_create_lattice(L, nd, size, "D", defaultNc);
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

/// QDP function profiling.
//  @function profile
//  @tparam[opt] integer level value to set (omit if just getting value)
//  @return previous value
static int
qopqdp_profile(lua_State *L)
{
  BEGIN_ARGS;
  OPT_INT(v, 0);
  END_ARGS;
  int r = QDP_profcontrol(v);
  if(nargs==0) {
    QDP_profcontrol(r);
  }
  lua_pushinteger(L, r);
  return 1;
}

static int
qopqdp_verbosity(lua_State *L)
{
  BEGIN_ARGS;
  OPT_INT(v, 0);
  END_ARGS;
  int qop = QOP_verbose(v);
  int qdp = QDP_verbose(v);
  if(nargs==0) {
    QOP_verbose(qop);
    QDP_verbose(qdp);
  }
  lua_pushinteger(L, qop);
  return 1;
}

static int
qopqdp_blocksize(lua_State *L)
{
  BEGIN_ARGS;
  OPT_INT(new, 0);
  END_ARGS;
  int old = QDP_get_block_size();
  lua_pushinteger(L, old);
  if(new>0) QDP_set_block_size(new);
  return 1;
}

static int
qopqdp_readGroupSize(lua_State *L)
{
  BEGIN_ARGS;
  OPT_INT(new, 0);
  END_ARGS;
  int old = QDP_set_read_group_size(new);
  lua_pushinteger(L, old);
  if(new<=0) QDP_set_read_group_size(old);
  return 1;
}

static int
qopqdp_writeGroupSize(lua_State *L)
{
  BEGIN_ARGS;
  OPT_INT(new, 0);
  END_ARGS;
  int old = QDP_set_write_group_size(new);
  lua_pushinteger(L, old);
  if(new<=0) QDP_set_write_group_size(old);
  return 1;
}

#if 0
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
#endif

static int
qopqdp_seed(lua_State *L)
{
  lattice_t *l = qopqdp_get_default_lattice(L);
  BEGIN_ARGS;
  GET_INT(seed);
  OPT_INT(uniform, -1);
  OPT_SUBSET(sub, l, QDP_all_L(l->qlat));
  END_ARGS;
  if(l->rs==NULL) {
    l->rs = QDP_create_S_L(l->qlat);
  }
  qhmc_qopqdp_seed_func(l->rs, seed, uniform, sub);
  qopqdp_srs = l->rs;
  QLA_Int i = 987654321 + QDP_this_node;
  QLA_S_eq_seed_i_I(&qopqdp_nrs, seed, &i);
  return 0;
}

static int
qopqdp_random(lua_State *L)
{
  BEGIN_ARGS;
  END_ARGS;
  double r = 0;
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
  OPT_INT(nc, 0);
  OPT_LATTICE(lat, NULL);
  END_ARGS;
  qopqdp_asqtad_create(L, nc, lat);
  return 1;
}

static int
qopqdp_hisq(lua_State *L)
{
  BEGIN_ARGS;
  OPT_INT(nc, 0);
  OPT_LATTICE(lat, NULL);
  END_ARGS;
  qopqdp_hisq_create(L, nc, lat);
  return 1;
}

static int
qopqdp_wilson(lua_State *L)
{
  BEGIN_ARGS;
  OPT_INT(nc, 0);
  OPT_LATTICE(lat, NULL);
  END_ARGS;
  qopqdp_wilson_create(L, nc, lat);
  return 1;
}

static int
qopqdp_dw(lua_State *L)
{
  BEGIN_ARGS;
  OPT_INT(nc, 0);
  OPT_LATTICE(lat, NULL);
  END_ARGS;
  qopqdp_dw_create(L, nc, lat);
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
// 3: (opt) lattice
static int
qopqdp_writer(lua_State *L)
{
  BEGIN_ARGS;
  GET_STRING(fn);
  GET_STRING(md);
  OPT_LATTICE(lat, NULL);
  END_ARGS;
  qopqdp_writer_create(L, fn, md, lat);
  return 1;
}

static int
qopqdp_getFileLattice(lua_State *L)
{
  BEGIN_ARGS;
  GET_STRING(fn);
  END_ARGS;
  QIO_Layout ql;
  ql.latdim = 0;
  ql.latsize = NULL;
  ql.this_node = QMP_get_node_number();
  ql.number_of_nodes = QMP_get_number_of_nodes();
  QIO_String *qs = QIO_string_create();
  QIO_Reader *qr = QIO_open_read(qs, fn, &ql, NULL, NULL);
  int nd = QIO_get_reader_latdim(qr);
  //printf0("lattice ndim = %i\n", *len);
  int *ls = QIO_get_reader_latsize(qr);
  //printf0("lattice size =");
  //for(i=0; i<nd; i++) printf0(" %i", ls[i]);
  //printf0("\n");
  qhmc_push_int_array(L, nd, ls);
  QIO_close_read(qr);
  return 1;
}

static int
qopqdp_remapout(lua_State *L)
{
  BEGIN_ARGS;
  GET_STRING(s);
  END_ARGS;
  // FIXME: make safe for multiple ranks
  int fd = creat(s, 0666);
  fflush(stdout);
  dup2(fd, 1);
  fflush(stderr);
  dup2(fd, 2);
  close(fd);
  return 0;
}

static int
qopqdp_option(lua_State *L)
{
  BEGIN_ARGS;
  GET_STRING(opt);
  if(nargs==2) nextarg++;
  END_ARGS;
  // return old value
#define GETI(s,v) if(strcmp(opt,s)==0) lua_pushinteger(L, v);
  GETI("MILCGaussian", QLA_use_milc_gaussian);
#undef GETI
  if(nargs==2) { // set new value
#define SETI(s,v) if(strcmp(opt,s)==0) v = luaL_checkinteger(L, 2);
    SETI("MILCGaussian", QLA_use_milc_gaussian);
#undef SETI
  }
  return 1;
}

static int
qopqdp_groupNumber(lua_State *L)
{
  BEGIN_ARGS;
  GET_STRING(s0);
  END_ARGS;
  int g = 0;
  const char *s = s0;
  while(s[0]!='\0') {
    switch(s[0]) {
    case 'T': g += GROUP_T; s++; break;
    case 'S': g += GROUP_S; s++; break;
    case 'G': g += GROUP_GL; s+=2; break;
    case 'U': g += GROUP_U; s++; break;
    case 'H': g += GROUP_H; s++; break;
    case 'A': g += GROUP_AH; s+=2; break;
    default: qlerror(L, 1, "unknown group string %s\n", s0);
    }
  }
  lua_pushinteger(L, g);
  return 1;
}

static int
qopqdp_groupName(lua_State *L)
{
  BEGIN_ARGS;
  GET_INT(g);
  END_ARGS;
  char s0[5];
  char *s = s0;
  if(g & GROUP_T) { s[0] = 'T'; s++; }
  if(g & GROUP_S) { s[0] = 'S'; s++; }
  switch(g&GROUP_TYPE) {
  case GROUP_GL: s[0]='G'; s[1]='L'; s+=2; break;
  case GROUP_U: s[0]='U'; s++; break;
  case GROUP_H: s[0]='H'; s++; break;
  case GROUP_AH: s[0]='A'; s[1]='H'; s+=2; break;
  default: qlerror(L, 1, "unknown group number %i\n", g);
  }
  s[0] = '\0';
  lua_pushstring(L, s0);
  return 1;
}

static struct luaL_Reg qopqdp_reg[] = {
  { "master",         qopqdp_master },
  { "rank",           qopqdp_rank },
  { "dtime",          qopqdp_dtime },
  { "defaultNc",      qopqdp_defaultNc },
  { "lattice",        qopqdp_lattice },
  { "defaultLattice", qopqdp_defaultLattice },
  { "profile",        qopqdp_profile },
  { "verbosity",      qopqdp_verbosity },
  { "blocksize",      qopqdp_blocksize },
  { "readGroupSize",  qopqdp_readGroupSize },
  { "writeGroupSize", qopqdp_writeGroupSize },
  { "seed",           qopqdp_seed },
  { "random",         qopqdp_random },
  { "gauge",          qopqdp_gauge },
  { "force",          qopqdp_force },
  { "asqtad",         qopqdp_asqtad },
  { "hisq",           qopqdp_hisq },
  { "wilson",         qopqdp_wilson },
  { "dw",             qopqdp_dw },
  { "reader",         qopqdp_reader },
  { "writer",         qopqdp_writer },
  { "getFileLattice", qopqdp_getFileLattice },
  { "remapout",       qopqdp_remapout },
  { "option",         qopqdp_option },
  { "groupNumber",    qopqdp_groupNumber },
  { "groupName",      qopqdp_groupName },
  { NULL, NULL}
};

void
qhmc_open_qopqdp(lua_State *L)
{
  luaL_register(L, "qopqdp", qopqdp_reg);
  int jobnum = QMP_get_job_number();
  int numjobs = QMP_get_number_of_jobs();
  lua_pushinteger(L, jobnum);
  lua_setglobal(L, "jobnum");
  lua_pushinteger(L, numjobs);
  lua_setglobal(L, "numjobs");
  lua_getglobal(L, "qopqdp");
  lua_pop(L, 1);
  QDP_set_read_group_size(64);
  QDP_set_write_group_size(64);
  open_qopqdp_smear(L);
#ifdef _OPENMP
  QDP_set_block_size(1024*1024);
#endif
}
