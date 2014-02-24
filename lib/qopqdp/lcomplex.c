#include "qhmc_qopqdp_common.h"

#define FTYPE complex
#define FTYPEC COMPLEX
#define ARITH
#define PRECISE
//#define COLORED
// QDP and QLA field type abbreviation
#define A C
#define AL c
// QDP field type name
#define QDPT Complex

#define GET_COLOR_SPIN (void)0
#define QLAELEM(x) (x)

#define LOOP_FTYPE_ELEM
#define END_LOOP_FTYPE_ELEM

#define GET_QLA_UNIT(x) QLA_Complex x; QLA_c_eq_r(x, 1)

static void
pushqlatype(lua_State *L, QLA_Complex c)
{
  qhmc_complex_create(L, QLA_real(c), QLA_imag(c));
}

#include "field0.c"
