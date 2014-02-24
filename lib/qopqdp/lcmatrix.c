#include "qhmc_qopqdp_common.h"

#define FTYPE cmatrix
#define FTYPEC CMATRIX
#define ARITH
#define PRECISE
#define COLORED
// QDP and QLA field type abbreviation
#define A M
#define AL m
// QDP field type name
#define QDPT ColorMatrix

#define GET_COLOR_SPIN GET_INT(ic); GET_INT(jc)
#define QLAELEM(x) QLA_elem_M(x,ic,jc)

#define LOOP_FTYPE_ELEM					\
  for(int ic=0; ic<QLA_Nc; ic++) {			\
  for(int jc=0; jc<QLA_Nc; jc++) {			\

#define END_LOOP_FTYPE_ELEM } }

#define GET_QLA_UNIT(x)					\
  QLA_ColorMatrix x;					\
  for(int ic=0; ic<QLA_Nc; ic++) {			\
    for(int jc=0; jc<QLA_Nc; jc++) {			\
      if(ic==jc) { QLA_c_eq_r(QLA_elem_M(x,ic,jc),1); }	\
      else { QLA_c_eq_r(QLA_elem_M(x,ic,jc),0);	}	\
    }							\
  }

static void
pushqlatype(lua_State *L, int NC, QLA_ColorMatrix m)
{
  //qhmc_complex_create(L, QLA_real(c), QLA_imag(c));
}

#include "field0.c"
