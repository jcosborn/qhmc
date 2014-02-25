#include "qhmc_qopqdp_common.h"
#include <math.h>

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
qlamakegroup(QLA_Complex *x, int g)
{
  switch(g&GROUP_TYPE) {
  case GROUP_GL: break;
  case GROUP_U: {
    QLA_Real n = QLA_norm2_c(*x);
    if(n==0) { QLA_c_eq_r(*x, 1); }
    else {
      n = 1/sqrt(n);
      QLA_c_eq_r_times_c(*x, n, *x);
    }
  } break;
  case GROUP_H: QLA_c_eq_r(*x, QLA_real(*x)); break;
  case GROUP_AH: QLA_c_eq_r_plus_ir(*x, 0, QLA_imag(*x)); break;
  }
  if(g&GROUP_S) QLA_c_eq_r(*x, 1);
  if(g&GROUP_T) QLA_c_eq_r(*x, 0);
}

static void
pushqlatype(lua_State *L, QLA_Complex *c)
{
  qhmc_complex_create(L, QLA_real(*c), QLA_imag(*c));
}

#include "field0.c"
