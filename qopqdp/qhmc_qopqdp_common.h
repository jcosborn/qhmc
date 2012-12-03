#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
#define QDP_Precision 'D'
#define QOP_Precision 'D'
#include <qop.h>
#include <qop_qdp.h>
#include "qhmc_qopqdp.h"

#if 0
#define lua_objlen lua_rawlen
#define luaL_register(L,n,l) do { \
    if(n) { luaL_newlibtable(L,l); lua_pushvalue(L,-1); lua_setglobal(L,n); } \
    luaL_setfuncs(L,l,0); \
  } while(0)
#endif

#define printf0 if(QDP_this_node==0) printf
#define printerr(...) fprintf(stderr, __VA_ARGS__)
#define ABORT(code) QDP_abort(code)
#define TRACE printf0("%s %s %i\n", __FILE__, __func__, __LINE__);

#define get_table_len(L, idx, n) luaL_checktype(L, idx, LUA_TTABLE); *(n) = lua_objlen(L, idx)
#define qerror(...) // FIXME
#define qassert(x) {\
    if(!(x)) {\
      TRACE;\
      printf0("assert failed (%s)\n", #x);\
      QDP_abort(1);\
    }\
  }

#define tableLoopKeys(L, idx)		\
  qassert(lua_type(L,idx)==LUA_TTABLE); \
  lua_pushnil(L);			\
  while(lua_next(L,idx)) {		\
    lua_pushvalue(L, -2);		\
    const char *tableLoopKey = lua_tostring(L, -1);
#define tableLoopKeysIfKeySetInt(L, key, var)			\
    if(strcmp(tableLoopKey,key)==0) (var) = lua_tointeger(L,-2)
#define tableLoopKeysIfKeySetDouble(L, key, var)			\
    if(strcmp(tableLoopKey,key)==0) (var) = lua_tonumber(L,-2)
#define tableLoopKeysEnd(L) \
  /*else {					\
      /-* error *-/				\
    }*/ \
    lua_pop(L, 2); \
  }
#define tableGetField(L, idx, key) lua_getfield(L, idx, key)


extern QLA_RandomState qopqdp_nrs;
extern QDP_RandomState *qopqdp_srs;

void get_bool_array(lua_State *L, int idx, int n, int *a);
void get_int_array(lua_State *L, int idx, int n, int *a);
void push_int_array(lua_State *L, int n, int *a);
void get_double_array(lua_State *L, int idx, int n, double *a);
void push_double_array(lua_State *L, int n, double *a);
const char *tableGetString(lua_State *L, int idx, char *key);
QOP_evenodd_t qopqdp_check_evenodd(lua_State *L, int idx);
QDP_Subset qopqdp_check_subset(lua_State *L, int idx);
QDP_Subset *qhmcqdp_get_timeslices(void);

QLA_Real infnorm_M(QDP_ColorMatrix *m, QDP_Subset s);
int check_uniform_M(QDP_ColorMatrix *m, QDP_Subset s);


typedef struct {
  QDP_ColorMatrix **links;
  QDP_ColorMatrix **lie;
  int nd;
} gauge_t;

gauge_t *qopqdp_gauge_create(lua_State *L);
gauge_t *qopqdp_gauge_check(lua_State *L, int idx);
void qopqdp_gauge_array_check(lua_State *L, int idx, int n, gauge_t *g[n]);


typedef struct {
  double time;
  double flops;
  int nd;
  QDP_ColorMatrix *force[];
} force_t;

force_t *qopqdp_force_create(lua_State *L);
force_t *qopqdp_force_check(lua_State *L, int idx);
void qopqdp_force_array_check(lua_State *L, int idx, int n, force_t *g[n]);
void check_force(force_t *f, gauge_t *g, double (*act)(gauge_t *g, void*), void*);
void get_lie_force(force_t *f, gauge_t *g);
void get_local_force(QLA_ColorMatrix *f, QLA_ColorMatrix *g0,
		     void (*update_link)(QLA_ColorMatrix *g, QLA_ColorMatrix *g0,
					 QLA_ColorMatrix *dg, double eps),
		     double (*act)(QLA_ColorMatrix *g, void *),
		     void *args, double eps);


void open_qopqdp_smear(lua_State* L);


typedef struct {
  double time;
  double flops;
  int its;
  QOP_hisq_coeffs_t coeffs;
  QOP_FermionLinksHisq *fl;
  QOP_F_FermionLinksHisq *ffl;
  gauge_t *g;
  double u0;
  double f7lf;
  QDP_Complex **fatphase, **longphase;
} hisq_t;

hisq_t *qopqdp_hisq_create(lua_State *L);
hisq_t *qopqdp_hisq_check(lua_State *L, int idx);
void hisqInvert(QOP_info_t *info, QOP_FermionLinksAsqtad *fla,
		QOP_invert_arg_t *invarg, QOP_resid_arg_t *residarg[],
		QLA_Real *masses, int nm, QDP_ColorVector *prop[],
		QDP_ColorVector *source);

typedef struct {
  QDP_ColorVector *cv;
} squark_t;

squark_t *qopqdp_squark_create(lua_State *L);
squark_t *qopqdp_squark_check(lua_State *L, int idx);
void qopqdp_squark_array_check(lua_State *L, int idx, int n, squark_t *q[n]);


typedef struct {
  double time;
  double flops;
  int its;
  QOP_wilson_coeffs_t coeffs;
  QOP_FermionLinksWilson *fl;
  QOP_F_FermionLinksWilson *ffl;
  gauge_t *g;
} wilson_t;

wilson_t *qopqdp_wilson_create(lua_State *L);
wilson_t *qopqdp_wilson_check(lua_State *L, int idx);

typedef struct {
  QDP_DiracFermion *df;
} wquark_t;

wquark_t *qopqdp_wquark_create(lua_State *L);
wquark_t *qopqdp_wquark_check(lua_State *L, int idx);
void qopqdp_wquark_array_check(lua_State *L, int idx, int n, wquark_t *q[n]);


typedef struct {
  double time;
  double flops;
  int its;
  int ls;
  QOP_dw_coeffs_t coeffs;
  QOP_FermionLinksDW *fl;
  QOP_F_FermionLinksDW *ffl;
  gauge_t *g;
} dw_t;

dw_t *qopqdp_dw_create(lua_State *L);
dw_t *qopqdp_dw_check(lua_State *L, int idx);

typedef struct {
  QDP_DiracFermion **df;
  int ls;
} dwquark_t;

dwquark_t *qopqdp_dwquark_create(lua_State *L, int ls);
dwquark_t *qopqdp_dwquark_check(lua_State *L, int idx);
void qopqdp_dwquark_array_check(lua_State *L, int idx, int n, dwquark_t *q[n]);
