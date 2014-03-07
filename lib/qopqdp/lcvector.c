#include "qhmc_qopqdp_common.h"
#include <math.h>
#include <qla_d.h>

#define FTYPE cvector
#define FTYPEC CVECTOR
#define ARITH
#define PRECISE
#define COLORED
// QDP and QLA field type abbreviation
#define A V
#define AL v
// QDP field type name
#define QDPT ColorVector

#define NREAL (2*(NC))
#define GET_COLOR_SPIN GET_INT(ic)
#define QLAELEM(x) QLA_elem_V(x,ic)

#define LOOP_FTYPE_ELEM	for(int ic=0; ic<QLA_Nc; ic++) {

#define END_LOOP_FTYPE_ELEM }

#define GET_QLA_CONST2(x,zr,zi)						\
  QLA_ColorVector x;							\
  for(int ic=0; ic<QLA_Nc; ic++) {					\
    QLA_c_eq_r_plus_ir(QLA_elem_V(x,ic),zr,zi);				\
  }

#define GET_QLA_UNIT(x)	   GET_QLA_CONST2(x,1,0)
#define GET_QLA_CONST(x,z) GET_QLA_CONST2(x,(z)->r,(z)->i)

// treat vector as diagonal matrix
static void
qlamakegroup(int NC, QLA_ColorVector *x, int g)
{
  switch(g&GROUP_TYPE) {
  case GROUP_GL: break;
  case GROUP_U: {
    for(int ic=0; ic<QLA_Nc; ic++) {
      QLA_Real n = QLA_norm2_c(QLA_elem_V(*x,ic));
      if(n==0) { QLA_c_eq_r(QLA_elem_V(*x,ic), 1); }
      else {
	n = 1/sqrt(n);
	QLA_c_eq_r_times_c(QLA_elem_V(*x,ic), n, QLA_elem_V(*x,ic));
      }
    }
  } break;
  case GROUP_H: {
    for(int ic=0; ic<QLA_Nc; ic++) {
      QLA_c_eq_r(QLA_elem_V(*x,ic), QLA_real(QLA_elem_V(*x,ic)));
    }
  } break;
  case GROUP_AH: {
    for(int ic=0; ic<QLA_Nc; ic++) {
      QLA_c_eq_r_plus_ir(QLA_elem_V(*x,ic), 0, QLA_imag(QLA_elem_V(*x,ic)));
    }
  } break;
  }
  if(g&GROUP_S) {
    QLA_D_Complex d1, d2;
    QLA_ColorVector y;
    QLA_Complex c;
    QLA_c_eq_r(c, 1);
    for(int ic=0; ic<QLA_Nc; ic++) {
      QLA_Complex z;
      QLA_c_eq_c_times_c(z, c, QLA_elem_V(*x,ic));
      QLA_c_eq_c(c, z);
    }
    QLA_c_eq_r_plus_ir(d1, QLA_real(c), QLA_imag(c));
    QLA_D_C_eq_clog_C(&d2, &d1);
    QLA_c_eq_r_times_c(d1, -1./QLA_Nc, d2);
    QLA_D_C_eq_cexp_C(&d2, &d1);
    QLA_c_eq_r_plus_ir(c, QLA_real(d2), QLA_imag(d2));
    QLA_V_eq_V(&y, x);
    QLA_V_eq_C_times_V(x, &c, &y);
  }
  if(g&GROUP_T) {
    QLA_Complex c1, c2;
    QLA_c_eq_r(c1, 0);
    for(int ic=0; ic<QLA_Nc; ic++) {
      QLA_c_peq_c(c1, QLA_elem_V(*x,ic));
    }
    QLA_c_eq_r_times_c(c2, 1./QLA_Nc, c1);
    for(int ic=0; ic<QLA_Nc; ic++) {
      QLA_c_meq_c(QLA_elem_V(*x,ic), c2);
    }
  }
}

static void
pushqlatype(lua_State *L, int NC, QLA_ColorVector *m)
{
  //qhmc_complex_create(L, QLA_real(c), QLA_imag(c));
}

#include "field0.c"
