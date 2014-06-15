#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

typedef struct {
  double r,i;
} qhmc_complex_t;

void open_qhmc_complex(lua_State* L);
qhmc_complex_t *qhmc_complex_create(lua_State* L, double re, double im);
qhmc_complex_t *qhmc_complex_check(lua_State *L, int idx);
qhmc_complex_t *qhmc_opt_complex(lua_State *L, int *idx, int required, qhmc_complex_t *def);
void qhmc_complex_get_as(lua_State *L, int idx, qhmc_complex_t *c);
qhmc_complex_t *qhmc_opt_as_complex_ptr(lua_State *L, int *idx, int required, qhmc_complex_t *buf, qhmc_complex_t *def);
void qhmc_opt_as_complex_def_real(lua_State *L, int *idx, int required, qhmc_complex_t *buf, double def);
int qhmc_opt_as_complex_array_len(lua_State *L, int idx, int required, int def);
void qhmc_opt_as_complex_array(lua_State *L, int *idx, int required, int n, qhmc_complex_t *t, int dn, qhmc_complex_t *def);

#define GET_COMPLEX(v) qhmc_complex_t *v=qhmc_complex_check(L,nextarg); nextarg++
#define OPT_COMPLEX(v,d) qhmc_complex_t *v=qhmc_opt_complex(L,&nextarg,0,d)
#define GET_AS_COMPLEX(v) qhmc_complex_t v; qhmc_complex_get_as(L,nextarg,&v); nextarg++
#define OPT_AS_COMPLEX_PTR(v,d) qhmc_complex_t _t ## v, *v = qhmc_opt_as_complex_ptr(L,&nextarg,0,&_t ## v,d)
#define OPT_AS_COMPLEX_DEF_REAL(v,d) qhmc_complex_t v; qhmc_opt_as_complex_def_real(L,&nextarg,0,&v,d)
#define OPT_AS_COMPLEX_ARRAY(n,t,dn,dt) int n=qhmc_opt_as_complex_array_len(L,nextarg,0,dn); qhmc_complex_t t[n==0?1:abs(n)]; qhmc_opt_as_complex_array(L,&nextarg,0,n,t,dn,dt)
