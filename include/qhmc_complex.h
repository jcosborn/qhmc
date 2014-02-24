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
