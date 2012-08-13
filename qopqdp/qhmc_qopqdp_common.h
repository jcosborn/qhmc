#include <stdlib.h>
#include <stdio.h>
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

extern QLA_RandomState qopqdp_nrs;
extern QDP_RandomState *qopqdp_srs;

extern void get_int_array(lua_State *L, int idx, int n, int *a);
extern void push_int_array(lua_State *L, int n, int *a);
extern void get_double_array(lua_State *L, int idx, int n, double *a);
extern void push_double_array(lua_State *L, int n, double *a);
extern QOP_evenodd_t qopqdp_check_evenodd(lua_State *L, int idx);
extern QDP_Subset qopqdp_check_subset(lua_State *L, int idx);
extern QDP_Subset *qhmcqdp_get_timeslices(void);

extern QLA_Real infnorm_M(QDP_ColorMatrix *m, QDP_Subset s);


typedef struct {
  QDP_ColorMatrix **links;
  QDP_ColorMatrix **lie;
  int nd;
} gauge_t;

extern gauge_t *qopqdp_gauge_create(lua_State* L);
extern gauge_t *qopqdp_gauge_check(lua_State *L, int idx);


typedef struct {
  double time;
  double flops;
  int nd;
  QDP_ColorMatrix *force[];
} force_t;

extern force_t *qopqdp_force_create(lua_State* L);
extern force_t *qopqdp_force_check(lua_State *L, int idx);
extern void check_force(force_t *f, gauge_t *g, double (*act)(gauge_t *g, void*), void*);
extern void get_lie_force(force_t *f, gauge_t *g);
extern void get_local_force(QLA_ColorMatrix *f, QLA_ColorMatrix *g0,
			    void (*update_link)(QLA_ColorMatrix *g, QLA_ColorMatrix *g0,
						QLA_ColorMatrix *dg, double eps),
			    double (*act)(QLA_ColorMatrix *g, void *),
			    void *args, double eps);


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

extern hisq_t *qopqdp_hisq_create(lua_State* L);
extern hisq_t *qopqdp_hisq_check(lua_State *L, int idx);
extern void
hisqInvert(QOP_info_t *info, QOP_FermionLinksAsqtad *fla, QOP_invert_arg_t *invarg,
	   QOP_resid_arg_t *residarg[], QLA_Real *masses, int nm,
	   QDP_ColorVector *prop[],  QDP_ColorVector *source);


typedef struct {
  QDP_ColorVector *cv;
} squark_t;

extern squark_t *qopqdp_squark_create(lua_State* L);
extern squark_t *qopqdp_squark_check(lua_State *L, int idx);
extern void qopqdp_squark_array_check(lua_State *L, int idx, int n, squark_t *q[n]);


typedef struct {
  double time;
  double flops;
  int its;
  QOP_wilson_coeffs_t coeffs;
  QOP_FermionLinksWilson *fl;
  QOP_F_FermionLinksWilson *ffl;
  gauge_t *g;
} wilson_t;

extern wilson_t *qopqdp_wilson_create(lua_State* L);
extern wilson_t *qopqdp_wilson_check(lua_State *L, int idx);


typedef struct {
  QDP_DiracFermion *df;
} wquark_t;

extern wquark_t *qopqdp_wquark_create(lua_State* L);
extern wquark_t *qopqdp_wquark_check(lua_State *L, int idx);
extern void qopqdp_wquark_array_check(lua_State *L, int idx, int n, wquark_t *q[n]);
