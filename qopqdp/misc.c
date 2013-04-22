#include "qhmc_qopqdp_common.h"
#include <string.h>
#include <math.h>
#include <qla_dn.h>
#include <qla_complex.h>
#include <complex.h>

#if 0
#if 0
#include <sys/types.h>
#include <sys/sysctl.h>
typedef long long timebase_t;
#define timer() __rdtsc()
#else
#include <mach/mach_time.h>
typedef uint64_t timebase_t;
#define timer() mach_absolute_time()
double timebaseSeconds(timebase_t tb)
{
  static double s = 0;
  if(s==0) {
    mach_timebase_info_data_t info;
    mach_timebase_info(&info);
    s = 1e-9*((double)info.numer)/((double)info.denom);
  }
  return s*tb;
}
#endif
#endif

#define NC nc
static void
infnorm(NCPROT QLA_ColorMatrix(*m), int i, void *args)
{
  QLA_Real *p = (QLA_Real*)args;
  QLA_Real r;
  QLA_r_eq_norm2_M(&r, m);
  if(r>*p) *p = r;
}
#undef NC

QLA_Real
infnorm_M(QDP_ColorMatrix *m, QDP_Subset s)
{
  QLA_Real nrm = 0;
  QDP_M_eq_funcia(m, infnorm, &nrm, s);
  double max = nrm;
  QMP_max_double(&max);
  return sqrt(max);
}

#define NC nc
static int isfirst=1;
static void
setfirst(NCPROT QLA_ColorMatrix(*m), int i, void *args)
{
  QLA_ColorMatrix(*r) = (QLA_ColorMatrix(*)) args;
  if(isfirst) {
    isfirst = 0;
    memcpy(r, m, sizeof(QLA_ColorMatrix()));
  }
}
#undef NC

#define NC nc
static int badIndexPlus1=0;
static void
compare(NCPROT QLA_ColorMatrix(*m), int i, void *args)
{
  QLA_ColorMatrix(*r) = (QLA_ColorMatrix(*))args;
  if(badIndexPlus1==0) {
    if(memcmp(m, r, sizeof(QLA_ColorMatrix()))) badIndexPlus1 = i+1;
  }
}
#undef NC

int
check_uniform_M(QDP_ColorMatrix *m, QDP_Subset s)
{
#define NC QDP_get_nc(m)
  QLA_ColorMatrix(r);
  isfirst = 1;
  QDP_M_eq_funcia(m, setfirst, &r, s);
  QMP_broadcast(&r, sizeof(QLA_ColorMatrix()));
  badIndexPlus1 = 0;
  QDP_M_eq_funcia(m, compare, &r, s);
  return badIndexPlus1;
#undef NC
}

//timebase_t tb0, tb1, tb2;

#define NC nc
// solves A^+ X + X B^+ = C for X
void
sylsolve_site(NCPROT QLA_ColorMatrix(*x), QLA_ColorMatrix(*a),
	      QLA_ColorMatrix(*b), QLA_ColorMatrix(*c))
{
  //timebase_t tb3 = timer();
  int nc1 = QOP_Nc;
  int nc2 = nc1*nc1;
  //printf("%i  %i\n", nc1, nc2);
  QLA_N_ColorMatrix(nc2, A);
  QLA_N_ColorVector(nc2, (*B)) = (QLA_N_ColorVector(nc2, (*))) c;
  QLA_N_ColorVector(nc2, (*X)) = (QLA_N_ColorVector(nc2, (*))) x;
  for(int i1=0; i1<nc1; i1++) {
    for(int i2=0; i2<nc1; i2++) {
      int i3 = nc1*i1 + i2;
      for(int k=0; k<nc2; k++) {
	QLA_c_eq_r(QLA_elem_M(A,i3,k), 0);
      }
      for(int k=0; k<nc1; k++) {
	int j3 = nc1*k + i2;  // i1 i2 k i2
	//QLA_c_peq_c(QLA_elem_M(A,i3,j3), QLA_elem_M(*a,i1,k));
	QLA_c_peq_ca(QLA_elem_M(A,i3,j3), QLA_elem_M(*a,k,i1));
	int jt = nc1*i1 + k;  // i1 i2 i1 k
	//QLA_c_peq_c(QLA_elem_M(A,i3,jt), QLA_elem_M(*b,k,i2));
	QLA_c_peq_ca(QLA_elem_M(A,i3,jt), QLA_elem_M(*b,i2,k));
      }
    }
  }
  //timebase_t tb4 = timer();
  //QLA_N_ColorMatrix(nc2, Ai);
#if QOP_Precision == 'F'
  //QLA_FN_M_eq_inverse_M(nc2, &Ai, &A);
  //QLA_FN_V_eq_M_times_V(nc2, X, &Ai, B);
  QLA_FN_V_eq_M_inverse_V(nc2, X, &A, B);
#else
  //QLA_DN_M_eq_inverse_M(nc2, &Ai, &A);
  //QLA_DN_V_eq_M_times_V(nc2, X, &Ai, B);
  QLA_DN_V_eq_M_inverse_V(nc2, X, &A, B);
#endif
  //timebase_t tb5 = timer();
  //tb0 += tb5 - tb3;
  //tb1 += tb4 - tb3;
  //tb2 += tb5 - tb4;
}
#undef NC

void
sqrt_deriv(QDP_ColorMatrix *deriv, QDP_ColorMatrix *sqrtM,
	   QDP_ColorMatrix *M, QDP_ColorMatrix *chain, QDP_Subset sub)
{
  int i;
  QDP_loop_sites(i, sub, {
      sylsolve_site(NCARG
		    (QLA_ColorMatrix(*))QDP_site_ptr_readwrite_M(deriv,i),
		    (QLA_ColorMatrix(*))QDP_site_ptr_readonly_M(sqrtM,i),
		    (QLA_ColorMatrix(*))QDP_site_ptr_readonly_M(sqrtM,i),
		    (QLA_ColorMatrix(*))QDP_site_ptr_readonly_M(chain,i));
    });
}

void
projectU_deriv(QDP_ColorMatrix *deriv, QDP_ColorMatrix *proj,
	       QDP_ColorMatrix *mat, QDP_ColorMatrix *chain, QDP_Subset sub)
{
  //tb0 = 0;
  //tb1 = 0;
  //tb2 = 0;
  //timebase_t tb3 = timer();
  int i;
  QDP_loop_sites(i, sub, {
      // F = C z - z (Cd U + Ud C) z (dy/du)
      QLA_ColorMatrix *d = (QLA_ColorMatrix(*))QDP_site_ptr_readwrite_M(deriv,i);
      QLA_ColorMatrix *p = (QLA_ColorMatrix(*))QDP_site_ptr_readonly_M(proj,i);
      QLA_ColorMatrix *m = (QLA_ColorMatrix(*))QDP_site_ptr_readonly_M(mat,i);
      QLA_ColorMatrix *c = (QLA_ColorMatrix(*))QDP_site_ptr_readonly_M(chain,i);
      //QLA_ColorMatrix x;
      QLA_ColorMatrix y;
      QLA_ColorMatrix z;
      //QLA_ColorMatrix cz;
      QLA_ColorMatrix t1;
      QLA_ColorMatrix t2;
      //QLA_M_eq_Ma_times_M(&x, m, m);
      QLA_M_eq_Ma_times_M(&y, m, p);
      QLA_M_eq_inverse_M(&z, &y);
      //QLA_M_eq_sqrt_M(&t1, &y);
      //QLA_M_eq_M_times_M(&cz, c, &z);
      //QLA_M_eq_M(d, &cz);
      //QLA_M_eq_Ma_times_M(&t1, p, &cz);
      QLA_M_eq_M_times_M(d, c, &z);
      QLA_M_eq_Ma_times_M(&t1, p, d);
      sylsolve_site(&t2, &y, &y, &t1);
      QLA_M_eq_M(&t1, &t2);
      QLA_M_peq_Ma(&t1, &t2);
      QLA_M_eq_M_times_M(&t2, m, &t1);
      QLA_M_meq_M(d, &t2);
    });
  //printf("%i  %i  %i\n", (int)tb0, (int)tb1, (int)tb2);
  //tb3 = timer() - tb3;
  //printf("projU: %g\n", timebaseSeconds(tb3));
}

void
exp_deriv_site(QLA_ColorMatrix *deriv, QLA_ColorMatrix *expM, 
	       QLA_ColorMatrix *M, QLA_ColorMatrix *chain) 
{
  // special SU(2) case
#if QDP_Colors == 2
  QLA_Complex Tr, det;
  QLA_c_eq_c_times_c(det, QLA_elem_M(*M,0,0),QLA_elem_M(*M,1,1));
  QLA_c_meq_c_times_c(det, QLA_elem_M(*M,0,1), QLA_elem_M(*M,1,0));

  QLA_C_eq_trace_M(&Tr, M);
  QLA_Complex qs;
  //  QLA_Complex qt;
  QLA_c_eq_r_times_c(qs, 0.5, Tr); // s=TrA/2

  QLA_Complex qs2;
  QLA_c_eq_c_times_c(qs2, qs, qs); //s2 = s^2
  QLA_c_meq_c(qs2, det); //s2 = s^2 - detA

  double _Complex t = QLA_real(qs2) + QLA_imag(qs2) * _Complex_I;
  t = csqrt(t); // sqrt(s^2 - det A)
  //  QLA_c_eq_r_plus_ir(qt, creal(t), cimag(t)); // t = sqrt(s^2 - det A)

  double _Complex exps, cosht, sinht, sinht_t;
  double _Complex s = QLA_real(qs) + QLA_imag(qs) * _Complex_I;
  exps = cexp(s);

  if(creal(t) == 0 && cimag(t) == 0) {
    cosht = 1;
    sinht = 0;
    sinht_t = 1;
  } else {
    cosht = ccosh(t);
    sinht = csinh(t);
    sinht_t = sinht/t;
  }

  double _Complex f0, f1;
  f1 = exps * sinht_t;
  f0 = exps * cosht - s * f1;;

  //derivative of the coefficients
  double _Complex f0s, f1s, f1t, f0t2, f1t2, t2;
  t2 = t*t;

  f0s = f0 - f1;
  f1s = f1;
  if (cabs(t) > 0.05) {
    f1t = ((f0-f1) + s*f1)/t;
    f1t2 = f1t/t;
  }
  else { //when |t| <= 0.05, the truncation error is O(10^{-14})
    f1t = exps * t/3 * (1+t2/10*(1+t2/28));
    f1t2 = exps / 3 * (1+t2/10*(1+t2/28));
  }
  f0t2 = f1 - s * f1t2;

  //QLA_Complex qf0;
  QLA_Complex qf1;
  //QLA_c_eq_r_plus_ir(qf0, creal(f0), cimag(f0));
  QLA_c_eq_r_plus_ir(qf1, creal(f1), cimag(f1));

  QLA_ColorMatrix B, AB;
  QLA_M_eq_M(&B, chain);
  QLA_M_eq_M_times_M(&AB, M, &B);
  QLA_M_eq_c_times_M(deriv, &qf1, &B); //f1 * B

  QLA_Complex trB, trAB;
  QLA_C_eq_trace_M(&trB,  &B);
  QLA_C_eq_trace_M(&trAB, &AB);

  double _Complex ctrB  = QLA_real(trB) + _Complex_I * QLA_imag(trB);
  double _Complex ctrAB = QLA_real(trAB) + _Complex_I * QLA_imag(trAB);
  double _Complex coeff;
  coeff  = (f0s - f0t2 * s) * ctrB;
  coeff += (f1s - f1t2 * s) * ctrAB;
  coeff *= 0.5;

  QLA_Complex qc;
  QLA_c_eq_r_plus_ir(qc, creal(coeff),cimag(coeff));

  QLA_M_peq_c(deriv, &qc); // f1 * B + () 
  
  coeff = 0.5 * (f0t2 * ctrB + f1t2 * ctrAB);
  QLA_c_eq_r_plus_ir(qc, creal(coeff),cimag(coeff));

  QLA_M_peq_c_times_M(deriv, &qc, M); // final result
  //  QLA_M_peq_Ma(deriv, deriv);

#else
  printerr("Not implemented\n");
  ABORT(1);
#endif
  
}

void
exp_deriv(QDP_ColorMatrix *deriv, QDP_ColorMatrix *expM,
	  QDP_ColorMatrix *M, QDP_ColorMatrix *chain, QDP_Subset sub)
{
  int i;
  QDP_loop_sites(i, sub, {
      exp_deriv_site(
		     (QLA_ColorMatrix(*))QDP_site_ptr_readwrite_M(deriv,i),
		     (QLA_ColorMatrix(*))QDP_site_ptr_readonly_M(expM,i),
		     (QLA_ColorMatrix(*))QDP_site_ptr_readonly_M(M,i),
		     (QLA_ColorMatrix(*))QDP_site_ptr_readonly_M(chain,i));
    });
}
