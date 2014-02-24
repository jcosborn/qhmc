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

#define GET_COLOR_SPIN (void)0
#define QLAELEM(x) (x)

#define LOOP_FTYPE_ELEM
#define END_LOOP_FTYPE_ELEM

#define GET_QLA_UNIT(x) QLA_Real x = 1

static void
pushqlatype(lua_State *L, QLA_Real r)
{
  lua_pushnumber(L, r);
}

#include "field0.c"
