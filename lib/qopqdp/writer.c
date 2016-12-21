#include "qhmc_qopqdp_common.h"

static char *mtname = "qopqdp.writer";

writer_t *
qopqdp_writer_check(lua_State *L, int idx)
{
  luaL_checkudata(L, idx, mtname);
  writer_t *w = lua_touserdata(L, idx);
#if 1
  int hasmt = lua_getmetatable(L, idx);
  qassert(hasmt==1);
  luaL_getmetatable(L, mtname);
  int eq = lua_equal(L, -1, -2);
  qassert(eq==1);
  lua_pop(L, 2);
#endif
  return w;
}

static int
qopqdp_writer_close(lua_State *L)
{
  BEGIN_ARGS;
  GET_WRITER(w);
  END_ARGS;
  if(w->open) {
    volatile int *p = &(w->pending);
    while(*p>0);
    QDP_close_write(w->qw);
    w->open = 0;
  }
  return 0;
}

#if 0
// 1: writer
// 2: field or table of fields of same type
// 3: metadata string
// 4: (need to add) optional options table
static int
qopqdp_writer_write(lua_State *L)
{
  int nargs = lua_gettop(L);
  qassert(nargs==3);
  writer_t *w = qopqdp_writer_check(L, 1);
  int nfields = 1;
  int istable = (lua_type(L,2)==LUA_TTABLE);
  if(istable) get_table_len(L, 2, &nfields);
  const char *mds = luaL_checkstring(L, 3);
  char prec = 'F';
  //QDP_set_write_group_size(8);
  //printf0("saving lattice file %s\n", fn);
  double dt = -QDP_time();
  QDP_String *md = QDP_string_create();
  QDP_string_set(md, (char *)mds);
  if(prec==QDP_Precision) {
    if(istable) {
      ldfermion_t *wq[nfields];
      qopqdp_ldfermion_array_check(L, 2, nfields, wq);
      QDP_DiracFermion *df[nfields];
      for(int i=0; i<nfields; i++) df[i] = wq[i]->df;
      QDP_vwrite_D(w->qw, md, df, nfields);
    } else {
      ldfermion_t *wq = qopqdp_ldfermion_check(L, 2);
      QDP_write_D(w->qw, md, wq->df);
    }
  } else {
#if QDP_Precision == 'F'
#define QDPO(x) QDP_D_ ## x
#define QDPOP(x) QDP_DF_ ## x
#else
#define QDPO(x) QDP_F_ ## x
#define QDPOP(x) QDP_FD_ ## x
#endif
    if(istable) {
#define NC QDP_get_nc(wq[0]->df)
      ldfermion_t *wq[nfields];
      qopqdp_ldfermion_array_check(L, 2, nfields, wq);
      QDPO(DiracFermion) *df[nfields];
      for(int i=0; i<nfields; i++) {
	df[i] = QDPO(create_D)();
	QDPOP(D_eq_D)(df[i], wq[i]->df, QDP_all);
      }
      QDPO(vwrite_D)(w->qw, md, df, nfields);
      for(int i=0; i<nfields; i++) {
	QDPO(destroy_D)(df[i]);
      }
#undef NC
    } else {
#define NC QDP_get_nc(wq->df)
      ldfermion_t *wq = qopqdp_ldfermion_check(L, 2);
      QDPO(DiracFermion) *df = QDPO(create_D)();
      QDPOP(D_eq_D)(df, wq->df, QDP_all);
      QDPO(write_D)(w->qw, md, df);
      QDPO(destroy_D)(df);
    }
#undef NC
#undef QDPO
#undef QDPOP
  }

  QDP_string_destroy(md);
  dt += QDP_time();
  printf0(" wrote in %g seconds\n", dt);
  return 0;
}
#endif

// 1: writer
// 2: table of Nc*Ns dfermion fields
// 3: metadata string
// 4: (need to add) optional options table
static int
qopqdp_writer_prop(lua_State *L)
{
  BEGIN_ARGS;
  GET_WRITER(w);
  GET_AS_QOPQDP_DFERMION_ARRAY(nfields, wq);
  GET_STRING(mds);
  OPT_STRING(prec, "F");
  END_ARGS;
  QDP_Lattice *qlat = w->lat->qlat;
  //QDP_set_write_group_size(8);
  //printf0("saving lattice file %s\n", fn);
  double dt = -QDP_time();
  QDP_String *md = QDP_string_create();
  QDP_string_set(md, (char *)mds);
#define NC QDP_get_nc(wq[0]->field)
  if(prec[0]==QDP_Precision) {
    qassert(nfields==(QLA_Nc*QLA_Ns));
    QDP_DiracPropagator *dp = QDP_create_P_L(qlat);
    for(int i=0; i<nfields; i++) {
      int color = i/QLA_Ns;
      int spin = i%QLA_Ns;
      QDP_P_eq_diracvec_D(dp, wq[i]->field, color, spin, QDP_all_L(qlat));
    }
    QDP_write_P(w->qw, md, dp);
    QDP_destroy_P(dp);
  } else {
#if QDP_Precision == 'F'
#define QDPO(x) QDP_D_ ## x
#define QDPOP(x) QDP_DF_ ## x
#else
#define QDPO(x) QDP_F_ ## x
#define QDPOP(x) QDP_FD_ ## x
#endif
    qassert(nfields==(QLA_Nc*QLA_Ns));
    QDPO(DiracPropagator) *dp = QDPO(create_P_L)(qlat);
    QDPO(DiracFermion) *df = QDPO(create_D_L)(qlat);
    for(int i=0; i<nfields; i++) {
      int color = i/QLA_Ns;
      int spin = i%QLA_Ns;
      QDPOP(D_eq_D)(df, wq[i]->field, QDP_all_L(qlat));
      QDPO(P_eq_diracvec_D)(dp, df, color, spin, QDP_all_L(qlat));
    }
    QDPO(write_P)(w->qw, md, dp);
    QDPO(destroy_D)(df);
    QDPO(destroy_P)(dp);
#undef QDPO
#undef QDPOP
  }
#undef NC

  QDP_string_destroy(md);
  dt += QDP_time();
  printf0(" saved in %g seconds\n", dt);
  return 0;
}

static struct luaL_Reg writer_reg[] = {
  { "__gc",    qopqdp_writer_close },
  { "close",   qopqdp_writer_close },
  //{ "write",   qopqdp_writer_write },
  { "prop",    qopqdp_writer_prop },
  { NULL, NULL}
};

writer_t *
qopqdp_writer_create(lua_State *L, const char *fn, const char *mds,
		     lattice_t *lat)
{
  if(lat==NULL) lat = qopqdp_get_default_lattice(L);
  writer_t *w = lua_newuserdata(L, sizeof(writer_t));
  QDP_String *md = QDP_string_create();
  QDP_string_set(md, (char *)mds);
  w->qw = QDP_open_write_L(lat->qlat, md, (char*)fn, QDP_SINGLEFILE);
  w->lat = lat;
  w->open = 1;
  w->pending = 0;
  QDP_string_destroy(md);
  if(luaL_newmetatable(L, mtname)) {
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, writer_reg);
  }
  lua_setmetatable(L, -2);
  return w;
}
