#include "qhmc_qopqdp_common.h"
#include <math.h>
#include <qla_d.h>

#define FTYPE dfermion
#define FTYPEC DFERMION
#define ISVECTOR
#define ARITH
#define PRECISE
#define COLORED
#define SPIN
// QDP and QLA field type abbreviation
#define A D
#define AL d
// QDP field type name
#define QDPT DiracFermion

#define NREAL ((NC)*QLA_Ns*2)
#define GET_COLOR_SPIN(i) GET_INT(S2(i,c)); GET_INT(S2(i,s))
#define OPT_COLOR_SPIN(i) OPT_INT(S2(i,c),-1); OPT_INT(S2(i,s),-1)
#define OPT_COLOR_SPIN_ZERO(i) OPT_INT(S2(i,c),0); OPT_INT(S2(i,s),0)
#define IS_SET_COLOR_SPIN(i) (S2(i,c)>=0 && S2(i,s)>=0)
#define CHECK_VALID_COLOR_SPIN(i) checkColorSpin(L,S2(i,c),S2(i,s),QLA_Nc)
#define QLAELEM(x,i) QLA_elem_D(x,S2(i,c),S2(i,s))
#define QLAELEMEQC(x,i,z) QLA_c_eq_r_plus_ir(QLAELEM(x,i),(z).r,(z).i)

#define LOOP_FTYPE_ELEM	\
  for(int ic=0; ic<QLA_Nc; ic++) { \
  for(int is=0; is<QLA_Ns; is++) {

#define END_LOOP_FTYPE_ELEM } }

#define GET_QLA_CONST2(x,zr,zi)						\
  QLA_DiracFermion(x);							\
  for(int ic=0; ic<QLA_Nc; ic++) {					\
    for(int is=0; is<QLA_Ns; is++) {					\
      QLA_c_eq_r_plus_ir(QLA_elem_D(x,ic,is),zr,zi);			\
    }									\
  }

#define GET_QLA_UNIT(x)	   GET_QLA_CONST2(x,1,0)
#define GET_QLA_CONST(x,z) GET_QLA_CONST2(x,(z)->r,(z)->i)

#define GET_QLA_CONST_ARRAY(x,z,n)					\
  qassert(n==(NC)*QLA_Ns);						\
  QLA_DiracFermion(x);							\
  for(int ic=0; ic<QLA_Nc; ic++) {					\
    for(int is=0; is<QLA_Ns; is++) {					\
      qhmc_complex_t *_pz = &(z)[ic*QLA_Ns+is];				\
      QLA_c_eq_r_plus_ir(QLA_elem_D(x,ic,is),_pz->r,_pz->i);		\
    }									\
  }

static void
checkColorSpin(lua_State *L, int ic, int is, int nc)
{
  if(ic<0 || ic>=nc) {
    qlerror(L, 1, "dirac fermion color index (%i) out of range 0..%i\n",
	    ic, nc-1);
  }
  if(is<0 || is>=4) {
    qlerror(L, 1, "dirac fermion spin index (%i) out of range 0..%i\n",
	    is, 3);
  }
}

static void
SPUR(QDP_Complex *d, QDP_DiracFermion *f, QDP_Subset sub)
{
#define NC QDP_get_nc(f)
  int s;
  QDP_loop_sites(s, sub, {
      QLA_Complex *ds = QDP_site_ptr_readwrite_C(d,s);
      QLA_DiracFermion(*fs) = QDP_site_ptr_readonly_D(f,s);
      QLA_Complex z;
      QLA_c_eq_r(z, 0);
      for(int ic=0; ic<QLA_Nc; ic++) {
	for(int is=0; is<QLA_Ns; is++) {
	  QLA_c_peq_c(z, QLA_elem_D(*fs,ic,is));
	}
      }
      QLA_c_eq_c(*ds, z);
    });
#undef NC
}

static void
DET(QDP_Complex *d, QDP_DiracFermion *f, QDP_Subset sub)
{
#define NC QDP_get_nc(f)
  int s;
  QDP_loop_sites(s, sub, {
      QLA_Complex *ds = QDP_site_ptr_readwrite_C(d,s);
      QLA_DiracFermion(*fs) = QDP_site_ptr_readonly_D(f,s);
      QLA_Complex z;
      QLA_Complex z2;
      QLA_c_eq_r(z, 1);
      for(int ic=0; ic<QLA_Nc; ic++) {
	for(int is=0; is<QLA_Ns; is++) {
	  QLA_c_eq_c_times_c(z2, z, QLA_elem_D(*fs,ic,is));
	  QLA_c_eq_c(z, z2);
	}
      }
      QLA_c_eq_c(*ds, z);
    });
#undef NC
}

// treat vector as diagonal matrix
static void
qlamakegroup(int NC, QLA_DiracFermion(*x), int g)
{
  switch(g&GROUP_TYPE) {
  case GROUP_GL: break;
  case GROUP_U: {
    for(int ic=0; ic<QLA_Nc; ic++) {
      for(int is=0; is<QLA_Ns; is++) {
	QLA_Real n = QLA_norm2_c(QLA_elem_D(*x,ic,is));
	if(n==0) { QLA_c_eq_r(QLA_elem_D(*x,ic,is), 1); }
	else {
	  n = 1/sqrt(n);
	  QLA_c_eq_r_times_c(QLA_elem_D(*x,ic,is), n, QLA_elem_D(*x,ic,is));
	}
      }
    }
  } break;
  case GROUP_H: {
    for(int ic=0; ic<QLA_Nc; ic++) {
      for(int is=0; is<QLA_Ns; is++) {
	QLA_c_eq_r(QLA_elem_D(*x,ic,is), QLA_real(QLA_elem_D(*x,ic,is)));
      }
    }
  } break;
  case GROUP_AH: {
    for(int ic=0; ic<QLA_Nc; ic++) {
      for(int is=0; is<QLA_Ns; is++) {
	QLA_c_eq_r_plus_ir(QLA_elem_D(*x,ic,is), 0, QLA_imag(QLA_elem_D(*x,ic,is)));
      }
    }
  } break;
  }
  if(g&GROUP_S) {
    QLA_D_Complex d1, d2;
    QLA_DiracFermion(y);
    QLA_Complex c;
    QLA_c_eq_r(c, 1);
    for(int ic=0; ic<QLA_Nc; ic++) {
      for(int is=0; is<QLA_Ns; is++) {
	QLA_Complex z;
	QLA_c_eq_c_times_c(z, c, QLA_elem_D(*x,ic,is));
	QLA_c_eq_c(c, z);
      }
    }
    QLA_c_eq_r_plus_ir(d1, QLA_real(c), QLA_imag(c));
    QLA_D_C_eq_clog_C(&d2, &d1);
    QLA_c_eq_r_times_c(d1, -1./(QLA_Nc*QLA_Ns), d2);
    QLA_D_C_eq_cexp_C(&d2, &d1);
    QLA_c_eq_r_plus_ir(c, QLA_real(d2), QLA_imag(d2));
    QLA_D_eq_D(&y, x);
    QLA_D_eq_C_times_D(x, &c, &y);
  }
  if(g&GROUP_T) {
    QLA_Complex c1, c2;
    QLA_c_eq_r(c1, 0);
    for(int ic=0; ic<QLA_Nc; ic++) {
      for(int is=0; is<QLA_Ns; is++) {
	QLA_c_peq_c(c1, QLA_elem_D(*x,ic,is));
      }
    }
    QLA_c_eq_r_times_c(c2, 1./(QLA_Nc*QLA_Ns), c1);
    for(int ic=0; ic<QLA_Nc; ic++) {
      for(int is=0; is<QLA_Ns; is++) {
	QLA_c_meq_c(QLA_elem_D(*x,ic,is), c2);
      }
    }
  }
}

static void
pushqlatype(lua_State *L, int NC, QLA_DiracFermion(*m))
{
  //qhmc_complex_create(L, QLA_real(c), QLA_imag(c));
}

#include "field0.c"
