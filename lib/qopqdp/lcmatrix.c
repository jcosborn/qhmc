#include "qhmc_qopqdp_common.h"
#include <qla_d.h>

#define FTYPE cmatrix
#define FTYPEC CMATRIX
#define ARITH
#define PRECISE
#define COLORED
#define COLORED2
// QDP and QLA field type abbreviation
#define A M
#define AL m
// QDP field type name
#define QDPT ColorMatrix

#define NREAL (2*(NC)*(NC))
#define GET_COLOR_SPIN(i) GET_INT(S2(i,c1)); GET_INT(S2(i,c2))
#define OPT_COLOR_SPIN(i) OPT_INT(S2(i,c1),-1); OPT_INT(S2(i,c2),-1)
#define OPT_COLOR_SPIN_ZERO(i) OPT_INT(S2(i,c1),0); OPT_INT(S2(i,c2),0)
#define IS_SET_COLOR_SPIN(i) (S2(i,c1)>=0 && S2(i,c2)>=0)
#define CHECK_VALID_COLOR_SPIN(i) checkColorSpin(L,S2(i,c1),S2(i,c2),QLA_Nc)
#define QLAELEM(x,i) QLA_elem_M(x,S2(i,c1),S2(i,c2))
#define QLAELEMEQC(x,i,z) QLA_c_eq_r_plus_ir(QLAELEM(x,i),(z).r,(z).i)
#define SPUR(d,f,s) QDP_C_eq_trace_M(d,f,s)
#define DET(d,f,s) QDP_C_eq_det_M(d,f,s)

#define LOOP_FTYPE_ELEM					\
  for(int ic1=0; ic1<QLA_Nc; ic1++) {			\
  for(int ic2=0; ic2<QLA_Nc; ic2++) {			\

#define END_LOOP_FTYPE_ELEM } }

#define GET_QLA_CONST2(x,zr,zi)						\
  QLA_ColorMatrix(x);							\
  for(int ic1=0; ic1<QLA_Nc; ic1++) {					\
    for(int ic2=0; ic2<QLA_Nc; ic2++) {					\
      if(ic1==ic2) { QLA_c_eq_r_plus_ir(QLA_elem_M(x,ic1,ic2),zr,zi); }	\
      else { QLA_c_eq_r(QLA_elem_M(x,ic1,ic2),0);	}		\
    }									\
  }

#define GET_QLA_UNIT(x)	   GET_QLA_CONST2(x,1,0)
#define GET_QLA_CONST(x,z) GET_QLA_CONST2(x,(z)->r,(z)->i)

#define GET_QLA_CONST_ARRAY(x,z,n)					\
  qassert(n==(NC)*(NC));						\
  QLA_ColorMatrix(x);							\
  for(int ic1=0; ic1<QLA_Nc; ic1++) {					\
    for(int ic2=0; ic2<QLA_Nc; ic2++) {					\
      qhmc_complex_t *_pz = &(z)[ic1*(NC)+ic2];				\
      QLA_c_eq_r_plus_ir(QLA_elem_M(x,ic1,ic2),_pz->r,_pz->i);		\
    }									\
  }

static void
checkColorSpin(lua_State *L, int ic1, int ic2, int nc)
{
  if(ic1<0 || ic1>=nc) {
    qlerror(L, 1, "first color matrix index (%i) out of range 0..%i\n",
	    ic1, nc-1);
  }
  if(ic2<0 || ic2>=nc) {
    qlerror(L, 1, "second color matrix index (%i) out of range 0..%i\n",
	    ic2, nc-1);
  }
}

static void
qlamakegroup(int NC, QLA_ColorMatrix(*x), int g)
{
  switch(g&GROUP_TYPE) {
  case GROUP_GL: break;
  case GROUP_U: {
    QLA_ColorMatrix(x2);
    QLA_ColorMatrix(n);
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
    QLA_ColorMatrix(y);
    QLA_M_eq_M(&y, x);
    QLA_M_peq_Ma(&y, x);
    QLA_Real half = 0.5;
    QLA_M_eq_r_times_M(x, &half, &y);
  } break;
  case GROUP_AH: {
    QLA_ColorMatrix(y);
    QLA_M_eq_M(&y, x);
    QLA_M_meq_Ma(&y, x);
    QLA_Real half = 0.5;
    QLA_M_eq_r_times_M(x, &half, &y);
  } break;
  }
  if(g&GROUP_T) {
    QLA_Complex c1, c2;
    QLA_C_eq_trace_M(&c1, x);
    QLA_c_eq_r_times_c(c2, 1./QLA_Nc, c1);
    QLA_M_meq_c(x, &c2);
  }
  if(g&GROUP_S) {
    QLA_D_Complex d1, d2;
    QLA_ColorMatrix(y);
    QLA_Complex c;
    QLA_C_eq_det_M(&c, x);
    QLA_c_eq_r_plus_ir(d1, QLA_real(c), QLA_imag(c));
    QLA_D_C_eq_clog_C(&d2, &d1);
    QLA_c_eq_r_times_c(d1, -1./QLA_Nc, d2);
    // FIX SH & SAH?
    //{
    //double w = (2*QHMC_PI)/QLA_Nc;
    //double im = QLA_imag(d1);
    //double im2 = im + round((QHMC_PI-im)/w)*w;
    QLA_D_C_eq_cexp_C(&d2, &d1);
    QLA_c_eq_r_plus_ir(c, QLA_real(d2), QLA_imag(d2));
    QLA_M_eq_M(&y, x);
    QLA_M_eq_C_times_M(x, &c, &y);
  }
}

static void
pushqlatype(lua_State *L, int NC, QLA_ColorMatrix(*m))
{
  int ns = NC*NC;
  qhmc_complex_t c[ns];
  for(int i=0; i<NC; i++) {
    for(int j=0; j<NC; j++) {
      QLA_Complex *z = &QLA_elem_M(*m,i,j);
      int k = i*NC + j;
      c[k].r = QLA_real(*z);
      c[k].i = QLA_imag(*z);
    }
  }
  qhmc_push_complex_array(L, ns, c);
}

#include "field0.c"
