#include "qhmc_qopqdp_common.h"
#include <qla_d.h>

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
qlamakegroup(int NC, QLA_ColorMatrix *x, int g)
{
  switch(g&GROUP_TYPE) {
  case GROUP_GL: break;
  case GROUP_U: {
    QLA_ColorMatrix x2, n;
    QLA_Complex d;
    QLA_M_eq_Ma_times_M(&x2, x, x);
    QLA_C_eq_det_M(&d, &x2);
    if(QLA_real(d)==0) {
      QLA_c_eq_r(d, 1);
      QLA_M_eq_c(x, &d);
    } else {
      QLA_M_eq_invsqrtPH_M(&n, &x2);
      QLA_M_eq_M_times_M(&x2, x, &n);
      QLA_M_eq_M(x, &x2);
    }
  } break;
  case GROUP_H: {
    QLA_ColorMatrix y;
    QLA_M_eq_M(&y, x);
    QLA_M_peq_Ma(&y, x);
    QLA_Real half = 0.5;
    QLA_M_eq_r_times_M(x, &half, &y);
  } break;
  case GROUP_AH: {
    QLA_ColorMatrix y;
    QLA_M_eq_M(&y, x);
    QLA_M_meq_Ma(&y, x);
    QLA_Real half = 0.5;
    QLA_M_eq_r_times_M(x, &half, &y);
  } break;
  }
  if(g&GROUP_S) {
    QLA_D_Complex d1, d2;
    QLA_ColorMatrix y;
    QLA_Complex c;
    QLA_C_eq_det_M(&c, x);
    QLA_c_eq_r_plus_ir(d1, QLA_real(c), QLA_imag(c));
    QLA_D_C_eq_clog_C(&d2, &d1);
    QLA_c_eq_r_times_c(d1, -1./QLA_Nc, d2);
    QLA_D_C_eq_cexp_C(&d2, &d1);
    QLA_c_eq_r_plus_ir(c, QLA_real(d2), QLA_imag(d2));
    QLA_M_eq_M(&y, x);
    QLA_M_eq_C_times_M(x, &c, &y);
  }
  if(g&GROUP_T) {
    QLA_Complex c1, c2;
    QLA_C_eq_trace_M(&c1, x);
    QLA_c_eq_r_times_c(c2, 1./QLA_Nc, c1);
    QLA_M_meq_c(x, &c2);
  }
}

static void
pushqlatype(lua_State *L, int NC, QLA_ColorMatrix *m)
{
  //qhmc_complex_create(L, QLA_real(c), QLA_imag(c));
}

#include "field0.c"
