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
  int nc1 = QLA_Nc;
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
#ifdef HAVENCN
#if QOP_Precision == 'F'
  //QLA_FN_M_eq_inverse_M(nc2, &Ai, &A);
  //QLA_FN_V_eq_M_times_V(nc2, X, &Ai, B);
  QLA_FN_V_eq_M_inverse_V(nc2, X, &A, B);
#else
  //QLA_DN_M_eq_inverse_M(nc2, &Ai, &A);
  //QLA_DN_V_eq_M_times_V(nc2, X, &Ai, B);
  QLA_DN_V_eq_M_inverse_V(nc2, X, &A, B);
#endif
#else
  qerror0(1,"can't use unitary projection without Nc=N libraries\n");
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
#define NC QDP_get_nc(deriv)
  int i;
  QDP_loop_sites(i, sub, {
      sylsolve_site(NCARG
		    (QLA_ColorMatrix(*))QDP_site_ptr_readwrite_M(deriv,i),
		    (QLA_ColorMatrix(*))QDP_site_ptr_readonly_M(sqrtM,i),
		    (QLA_ColorMatrix(*))QDP_site_ptr_readonly_M(sqrtM,i),
		    (QLA_ColorMatrix(*))QDP_site_ptr_readonly_M(chain,i));
    });
#undef NC
}

void
projectU_deriv(QDP_ColorMatrix *deriv, QDP_ColorMatrix *proj,
	       QDP_ColorMatrix *mat, QDP_ColorMatrix *chain, QDP_Subset sub)
{
#define NC QDP_get_nc(deriv)
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
      sylsolve_site(NCARG &t2, &y, &y, &t1);
      QLA_M_eq_M(&t1, &t2);
      QLA_M_peq_Ma(&t1, &t2);
      QLA_M_eq_M_times_M(&t2, m, &t1);
      QLA_M_meq_M(d, &t2);
    });
  //printf("%i  %i  %i\n", (int)tb0, (int)tb1, (int)tb2);
  //tb3 = timer() - tb3;
  //printf("projU: %g\n", timebaseSeconds(tb3));
#undef NC
}

#ifdef HAVENC3
#define NC 3
void 
traceless_herm_M_evalues(QLA_ColorMatrix *Q, double _Complex *u, double _Complex *w, 
			   double _Complex *q1, double _Complex *q2, double _Complex *q3) {
  
  QLA_Complex c0, c1;

  QLA_ColorMatrix Q2;
  QLA_M_eq_M_times_M(&Q2, Q, Q);

  QLA_C_eq_det_M    (&c0, Q);       // c0 = det(Q)
  QLA_C_eq_trace_M  (&c1, &Q2);     // c1 = tr(Q^2)
  
  double athird = 1.0/3.0;
  double _Complex cc0, cc1, cc0max;
  QLA_c99_eq_c(cc0, c0);
  QLA_c99_eq_c(cc1, c1);
  cc1 *= 0.5;

  cc0max = 2*csqrt(cc1 * athird)*(cc1 * athird);  //c_0^max = 2 * (c1/3)^{3/2}
  double _Complex theta;
  
  theta =cacos(cc0/cc0max);
  *u = csqrt(athird * cc1) * ccos(athird * theta);
  *w = csqrt(cc1) * csin(athird * theta);
  *q1 = 2 * *u;
  *q2 = -*u + *w;
  *q3 = -*u - *w;
}

static void
get_Bs(QLA_ColorMatrix *Q, QLA_ColorMatrix *Q2, QLA_ColorMatrix *B1,
       QLA_ColorMatrix *B2, double _Complex *f0, double _Complex *f1,
       double _Complex *f2)
{

  double _Complex u, w, q1, q2, q3;
  traceless_herm_M_evalues(Q, &u, &w, &q1, &q2, &q3);

  double _Complex e2iu, e_iu;
  e2iu = cexp(2 * _Complex_I * u);
  e_iu = cexp(-1.0 * _Complex_I * u);
  
  double _Complex u2 = u*u;
  double _Complex w2 = w*w;

  double _Complex zeta0w, zeta1w;
  if (cabs(w) > 0.05) {
    zeta0w = sin(w)/w;
    zeta1w = (cos(w)-zeta0w)/w2;
  }
  else {
    zeta0w = 1 - w2/6. * (1-w2/20. * (1 - w2/42.));
    zeta1w = -(1 - w2/10. * (1 - w2/28.*(1 - w2/54.)))/3.0;
  }
  
  double _Complex h0, h1, h2; 
  h0 = (u2 - w2) * e2iu + e_iu * ( 8 * u2 *cos(w) + 2*_Complex_I*u * (3*u2+w2)*zeta0w);
  h1 = 2*u*e2iu - e_iu * (2 * u * cos(w) - _Complex_I * (3*u2-w2)*zeta0w);
  h2 = e2iu - e_iu * ( cos(w) + 3*_Complex_I*u * zeta0w);

  double _Complex fac = 1.0/(9*u2-w2);
  *f0 = h0 * fac;
  *f1 = h1 * fac;
  *f2 = h2 * fac;

  double _Complex r01, r11, r21, r02, r12, r22, iu;
  double _Complex cosw = ccos(w);

  iu = _Complex_I * u;

  r01 = 2*(u + _Complex_I * (u2 - w2)) * e2iu 
	   + 2 * e_iu * ( 4*u*(2 - iu) * cosw + _Complex_I * (9 * u2 + w2 - iu * (3*u2 + w2))*zeta0w);

  r11 = 2*(1 + 2*iu) * e2iu 
    + e_iu * ( -2 * (1-iu) * cosw + _Complex_I * (6*u + _Complex_I * (w2 - 3*u2)) * zeta0w);

  r21 = 2 * _Complex_I * e2iu + _Complex_I * e_iu * (cosw - 3*(1-iu)*zeta0w);

  r02 = -2 * e2iu + 2 * iu * e_iu * (cosw + (1+4*iu) * zeta0w + 3 * u2 * zeta1w);

  r12 = -_Complex_I * e_iu * ( cosw + (1+2*iu) * zeta0w - 3*u2 * zeta1w);

  r22 = e_iu * (zeta0w - 3 * iu * zeta1w);

  double _Complex b10, b11, b12, b20, b21, b22;
  
  double _Complex fac1, fac2, fac3;
  
  double _Complex mult = 0.5 * fac * fac;
  fac1 = 2 * u;
  fac2 = 3*u2 - w2;
  fac3 = 2*(15*u2+w2);

  b10 = fac1 * r01 + fac2 * r02 - fac3 * (*f0);
  b10 *= mult;

  b11 = fac1 * r11 + fac2 * r12 - fac3 * (*f1);
  b11 *= mult;

  b12 = fac1 * r21 + fac2 * r22 - fac3 * (*f2);
  b12 *= mult;

  fac2 = 3*u;
  fac3 = 24*u;
  b20 = r01 - fac2 * r02 - fac3 * (*f0);
  b20 *= mult;

  b21 = r11 - fac2 * r12 - fac3 * (*f1);
  b21 *= mult;

  b22 = r21 - fac2 * r22 - fac3 * (*f2);
  b22 *= mult;

  QLA_Complex qb10, qb11, qb12, qb20, qb21, qb22;
  QLA_c_eq_c99(qb10, b10);
  QLA_c_eq_c99(qb11, b11);
  QLA_c_eq_c99(qb12, b12);
  QLA_c_eq_c99(qb20, b20);
  QLA_c_eq_c99(qb21, b21);
  QLA_c_eq_c99(qb22, b22);

  QLA_M_eq_c(B1, &qb10);
  QLA_M_peq_c_times_M(B1, &qb11, Q); 
  QLA_M_peq_c_times_M(B1, &qb12, Q2);
  
  QLA_M_eq_c(B2, &qb20);
  QLA_M_peq_c_times_M(B2, &qb21, Q); 
  QLA_M_peq_c_times_M(B2, &qb22, Q2);
}
#undef NC
#endif

#define NC nc
void
exp_deriv_site(NCPROT QLA_ColorMatrix *deriv, QLA_Real *r, 
	       QLA_ColorMatrix *M, QLA_ColorMatrix *chain)
{
  QLA_ColorMatrix tmp;
  QLA_ColorMatrix A;
  QLA_M_eq_r_times_M(&A, r, M);
  
#ifdef HAVENC3
  // special SU(3) case
  if(QLA_Nc==3) {
    QLA_Complex minus_i;
    QLA_c_eq_r_plus_ir(minus_i, 0, -1);

    QLA_ColorMatrix Q, Q2;
    QLA_M_eq_C_times_M(&Q, &minus_i, &A);
    QLA_M_eq_M_times_M(&Q2, &Q, &Q);

    double _Complex f0, f1, f2;
    QLA_ColorMatrix B1, B2;
    get_Bs(&Q, &Q2, &B1, &B2, &f0, &f1, &f2);

    QLA_Complex /*qf0,*/ qf1, qf2;
    //QLA_c_eq_c99(qf0, f0);
    QLA_c_eq_c99(qf1, f1);
    QLA_c_eq_c99(qf2, f2);
 
    // derivative  
    QLA_Complex trB1M, trB2M;
    QLA_ColorMatrix prod, mat;
    QLA_M_eq_Ma (&mat, chain);

    //tr(B_1 M)
    QLA_M_eq_M_times_M (&prod, &B1, &mat); //B_1 M
    QLA_C_eq_trace_M   (&trB1M, &prod);
  
    //tr(B_2 M);
    QLA_M_eq_M_times_M (&prod, &B2, &mat); //B_2 M
    QLA_C_eq_trace_M   (&trB2M, &prod);
  
    // deriv = Tr(B_1 M) Q
    QLA_M_eq_c_times_M (&tmp, &trB1M, &Q);
  
    // deriv += Tr(B_2 M) Q^2
    QLA_M_peq_c_times_M (&tmp, &trB2M, &Q2);
  
    // deriv += f1 M
    QLA_M_peq_c_times_M (&tmp, &qf1, &mat);
  
    // deriv += f2 Q M
    QLA_M_eq_M_times_M  (&prod, &Q, &mat); // Q M
    QLA_M_peq_c_times_M (&tmp, &qf2, &prod); 
  
    // deriv += f2 M Q
    QLA_M_eq_M_times_M  (&prod, &mat, &Q); // M Q
    QLA_M_peq_c_times_M (&tmp, &qf2, &prod);

    QLA_M_eq_c_times_M  (&tmp, &minus_i, &tmp);

    QLA_M_eq_Ma(deriv, &tmp);
  }
#endif
#ifdef HAVENC2
  // special SU(2) case
  if(QLA_Nc==2) {
    QLA_Complex Tr, det;
    QLA_c_eq_c_times_c(det, QLA_elem_M(A,0,0),QLA_elem_M(A,1,1));
    QLA_c_meq_c_times_c(det, QLA_elem_M(A,0,1), QLA_elem_M(A,1,0));

    QLA_C_eq_trace_M(&Tr, &A);

    if (QLA_real(Tr) == 0 && QLA_imag(Tr) == 0) {
      double _Complex t, t2;
      double _Complex cosht, sinht, sinht_t;

      t = csqrt(-QLA_real(det) - _Complex_I * QLA_imag(det)); 
      t2 = t * t;

      if(creal(t) == 0 && cimag(t) == 0) {
	cosht = 1;
	sinht = 0;
	sinht_t = 1;
      } else {
	cosht = ccosh(t);
	sinht = csinh(t);
	sinht_t = sinht / t;
      }

      double _Complex f0, f1;
      f0 = cosht;
      f1 = sinht_t;

      double _Complex f1t, f0t2, f1t2;
      if (cabs(t) > 0.05) {
	f1t = (f0-f1) / t;
	f1t2 = f1t / t;
      }
      else { //when |t| <= 0.05, the truncation error is O(10^{-14})
	f1t = t / 3.0 * (1+t2/10*(1+t2/28));
	f1t2 = 1.0 / 3.0 * (1+t2/10*(1+t2/28));
      }
      f0t2 = f1;

      QLA_Complex qf1;
      QLA_ColorMatrix B, AB;
      QLA_c_eq_r_plus_ir(qf1, creal(f1), cimag(f1));

      QLA_M_eq_Ma(&B, chain);
      QLA_M_eq_M_times_M(&AB, &A, &B);
      QLA_M_eq_c_times_M(&tmp, &qf1, &B); //f1 * B

      QLA_Complex trB, trAB;
      QLA_C_eq_trace_M(&trB,  &B);
      QLA_C_eq_trace_M(&trAB, &AB);

      double _Complex ctrB  = QLA_real(trB) + _Complex_I * QLA_imag(trB);
      double _Complex ctrAB = QLA_real(trAB) + _Complex_I * QLA_imag(trAB);
      double _Complex coeff;

      coeff = 0.5 * (f0t2 * ctrB + f1t2 * ctrAB);
      QLA_Complex qc;
      QLA_c_eq_r_plus_ir(qc, creal(coeff),cimag(coeff));
      QLA_M_peq_c_times_M(&tmp, &qc, &A);
    } else {
      QLA_Complex qs;
      QLA_c_eq_r_times_c(qs, 0.5, Tr); // s=TrA/2

      QLA_Complex qs2;
      QLA_c_eq_c_times_c(qs2, qs, qs); //s2 = s^2
      QLA_c_meq_c(qs2, det); //s2 = s^2 - detA

      double _Complex t = QLA_real(qs2) + QLA_imag(qs2) * _Complex_I;
      t = csqrt(t); // sqrt(s^2 - det A)

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
      double _Complex f0s, f1s, f0t, f1t, f0t2, f1t2, t2;
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
      f0t  = t*f1 - s*f1t;
      f0t2 = f1 - s*f1t2;

      QLA_Complex qf1;
      QLA_c_eq_r_plus_ir(qf1, creal(f1), cimag(f1));

      QLA_ColorMatrix B, AB;
      QLA_M_eq_Ma(&B, chain);
      QLA_M_eq_M_times_M(&AB, &A, &B);
      QLA_M_eq_c_times_M(&tmp, &qf1, &B); //f1 * B

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

      QLA_M_peq_c(&tmp, &qc); // f1 * B + () 

      coeff = 0.5 * (f0t2 * ctrB + f1t2 * ctrAB);

      QLA_c_eq_r_plus_ir(qc, creal(coeff),cimag(coeff));
      QLA_M_peq_c_times_M(&tmp, &qc, &A);
    }
    QLA_M_eq_Ma(deriv, &tmp);
  }
#endif
  if(QLA_Nc!=2 || QLA_Nc!=3) {
    qerror0(1,"Not implemented\n");
  }
}
#undef NC

void
exp_deriv(QDP_ColorMatrix *deriv, QLA_Real *r,
	  QDP_ColorMatrix *M, QDP_ColorMatrix *chain, QDP_Subset sub)
{
#define NC QDP_get_nc(deriv)
  int i;
  QDP_loop_sites(i, sub, {
      exp_deriv_site(NCARG
		     (QLA_ColorMatrix(*))QDP_site_ptr_readwrite_M(deriv,i),
		     r,
		     (QLA_ColorMatrix(*))QDP_site_ptr_readonly_M(M,i),
		     (QLA_ColorMatrix(*))QDP_site_ptr_readonly_M(chain,i));
    });
#undef NC
}
