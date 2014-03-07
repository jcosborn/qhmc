#include "qhmc_qopqdp_common.h"

#define FTYPE real
#define FTYPEC REAL
#define ISREAL
#define ARITH
#define PRECISE
//#define COLORED
// QDP and QLA field type abbreviation
#define A R
#define AL r
// QDP field type name
#define QDPT Real

#define NREAL 1
#define GET_COLOR_SPIN (void)0
#define QLAELEM(x) (x)

#define LOOP_FTYPE_ELEM
#define END_LOOP_FTYPE_ELEM

#define GET_QLA_UNIT(x) QLA_Real x = 1
#define GET_QLA_CONST(x,z) QLA_Real x = (z)->r

static void
qlamakegroup(QLA_Real *x, int g)
{
  switch(g&GROUP_TYPE) {
  case GROUP_GL: break;
  case GROUP_U: *x = (*x>=0) ? 1 : -1; break;
  case GROUP_H: break;
  case GROUP_AH: *x = 0; break;
  }
  if(g&GROUP_S) *x = 1;
  if(g&GROUP_T) *x = 0;
}

static void
pushqlatype(lua_State *L, QLA_Real *r)
{
  lua_pushnumber(L, *r);
}

#include "field0.c"
