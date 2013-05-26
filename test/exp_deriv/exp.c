#include <qhmc_qopqdp_common.h>
#include <complex.h>
#include <math.h>
#include <qla_complex.h>

void
printm(QLA_ColorMatrix *m)
{
  for(int i=0; i<QLA_Nc; i++) {
    for(int j=0; j<QLA_Nc; j++) {
      QLA_Complex *z = &QLA_elem_M(*m,i,j);
      printf("%10g%+10gi  ", QLA_real(*z), QLA_imag(*z));
    }
    printf("\n");
  }
  printf("\n");
}

void
printc(QLA_Complex *c) {
  printf("%f+i%f\n",QLA_real(*c), QLA_imag(*c));
}

void
printc99(double _Complex *c) {
  printf("%f+i%f\n",creal(*c), cimag(*c));
}

/*
void 
traceless_herm_M_evalues(QLA_ColorMatrix *Q, double *u, double *w, 
			   double *q1, double *q2, double *q3) {
  
  QLA_Complex c0, c1;

  QLA_ColorMatrix Q2;
  QLA_M_eq_M_times_M(&Q2, Q, Q);
  printf("Q^2 = \n"); printm(&Q2);

  QLA_C_eq_det_M    (&c0, Q);     // c0 = det(Q)
  QLA_C_eq_trace_M  (&c1, &Q2);     // c1 = tr(Q^2)
  

  double athird = 1.0/3.0;
  double cc0, cc1, cc0max;
  cc0 = QLA_real(c0);     
  cc1 = 0.5 * QLA_real(c1);
  cc0max = 2*sqrt(cc1 * athird)*(cc1 * athird);  //c_0^max = 2 * (c1/3)^{3/2}
  printf("c0 = %f\n", cc0);
  printf("c1 = %f\n", cc1);
  printf("c0_max = %f\n", cc0max);

  double theta;
  
  theta =acos(cc0/cc0max);
  *u = sqrt(athird * cc1) * cos(athird * theta);
  *w = sqrt(cc1) * sin(athird * theta);
  *q1 = 2 * *u;
  *q2 = -*u + *w;
  *q3 = -*u - *w;

  printf("u = %f, w = %f, q1 = %f, q2 = %f, q3 = %f\n", *u, *w, *q1, *q2, *q3);
}

void 
get_f_coeffs(QLA_ColorMatrix *Q, double _Complex *f0, double _Complex *f1, double _Complex *f2){
  double u, w, q1, q2, q3;
  traceless_herm_M_evalues(Q, &u, &w, &q1, &q2, &q3);
  printf("q1=\n"); printc99(&q1);

  double _Complex e2iu, e_iu;
  e2iu = cexp(2 * _Complex_I * u);
  e_iu = cexp(-1.0 * _Complex_I * u);
  
  double u2 = u*u;
  double w2 = w*w;

  double _Complex zeta0w;
  if (fabs(w) > 0.05) {
    zeta0w = sin(w)/w;
  }
  else {
    zeta0w = 1 - w2/6. * (1-w2/20. * (1 - w2/42.));
  }
  
  double _Complex h0, h1, h2;
  h0 = (u2 - w2) * e2iu + e_iu * ( 8 * u2 *cos(w) + 2*_Complex_I*u * (3*u2+w2)*zeta0w);
  h1 = 2*u*e2iu - e_iu * (2 * u * cos(w) - _Complex_I * (3*u2-w2)*zeta0w);
  h2 = e2iu - e_iu * ( cos(w) + 3*_Complex_I*u * zeta0w);

  double fac = 1.0/(9*u2-w2);
  *f0 = h0 * fac;
  *f1 = h1 * fac;
  *f2 = h2 * fac;
}

void
get_Bs(QLA_ColorMatrix *Q, QLA_ColorMatrix *Q2, QLA_ColorMatrix *B1, QLA_ColorMatrix *B2, double _Complex *f0, double _Complex *f1, double _Complex *f2) {
  double u, w, q1, q2, q3;
  traceless_herm_M_evalues(Q, &u, &w, &q1, &q2, &q3);
  printf("q1=\n"); printc99(&q1);

  double _Complex e2iu, e_iu;
  e2iu = cexp(2 * _Complex_I * u);
  e_iu = cexp(-1.0 * _Complex_I * u);
  
  double u2 = u*u;
  double w2 = w*w;

  double _Complex zeta0w, zeta1w;
  if (fabs(w) > 0.05) {
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

  double fac = 1.0/(9*u2-w2);
  *f0 = h0 * fac;
  *f1 = h1 * fac;
  *f2 = h2 * fac;

  double _Complex r01, r11, r21, r02, r12, r22, iu;
  double cosw = cos(w);

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
  
  double fac1, fac2, fac3;
  
  double mult = 0.5 * fac * fac;
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

*/
int
main(void) {
  QLA_ColorMatrix O, iQ, matI;
  QLA_M_eq_zero(&matI);

  for(int i=0; i<QLA_Nc; i++) {
    QLA_c_eq_r_plus_ir(QLA_elem_M(matI,i,i), 1.0, 0.0);
  }

  printm(&matI);

  //QLA_Complex tr;
  //QLA_Real half = 0.5;

  for(int i=0; i<QLA_Nc; i++) {
    for(int j=0; j<QLA_Nc; j++) {
      QLA_c_eq_r_plus_ir(QLA_elem_M(O,i,j), i+1, QLA_Nc*(j+1));
    }
    QLA_c_eq_r_plus_ir(QLA_elem_M(O,i,i), 2+1, 1);
  }
  printm(&O);

#if QDP_Colors == 3
  QLA_ColorMatrix A;
  QLA_M_eq_zero(&A);

  for ( int m = 0; m < QLA_Nc; m++) {
    for ( int n = 0; n < QLA_Nc; n++) {
      QLA_c_eq_r_plus_ir(QLA_elem_M(A, m, n), 3+m, 2-n);
    }
    QLA_c_eq_r_plus_ir(QLA_elem_M(A,m,m), 2+1, 1);

  }
  QLA_M_eq_antiherm_M(&A, &A);
  printm(&A);
  
  QLA_M_eq_zero(&A);

  QLA_c_eq_r_plus_ir(QLA_elem_M(A,0,0),0,-1);
  QLA_c_eq_r_plus_ir(QLA_elem_M(A,1,1),0,0.4);
  QLA_c_eq_r_plus_ir(QLA_elem_M(A,2,2),0,0.6);
  printm(&A);
  
  QLA_Complex minus_i;
  QLA_c_eq_r_plus_ir(minus_i, 0, -1);
  QLA_ColorMatrix Q, Q2, expiQ, qla_expA;

  QLA_M_eq_C_times_M(&Q, &minus_i, &matI);
  QLA_M_eq_M_times_M(&Q2, &Q, &Q);

  printf("Q=\n"); printm(&Q);

  double _Complex f0, f1, f2;
  QLA_ColorMatrix B1, B2;
  
  get_Bs(&Q, &Q2, &B1, &B2, &f0, &f1, &f2);
  printf("f0, f1, f2=\n");
  printc99(&f0);
  printc99(&f1);
  printc99(&f2);

  QLA_Complex qf0, qf1, qf2;

  QLA_c_eq_c99(qf0, f0);
  QLA_c_eq_c99(qf1, f1);
  QLA_c_eq_c99(qf2, f2);
    
  QLA_M_eq_c(&expiQ, &qf0);
  QLA_M_peq_c_times_M(&expiQ, &qf1, &Q);
  QLA_M_peq_c_times_M(&expiQ, &qf2, &Q2);

  QLA_M_eq_exp_M(&qla_expA, &matI);

  printf("my expiQ = \n");
  printm(&expiQ);
  
  printf("qla expA = \n");
  printm(&qla_expA);

  // derivative
   
  QLA_Complex trB1M, trB2M;
  QLA_ColorMatrix prod, deriv;

  //tr(B_1 M)
  QLA_M_eq_M_times_M (&prod, &B1, &matI); //B_1 M
  QLA_C_eq_trace_M   (&trB1M, &prod);
  
  //tr(B_2 M);
  QLA_M_eq_M_times_M (&prod, &B2, &matI); //B_2 M
  QLA_C_eq_trace_M   (&trB2M, &prod);
  
  // deriv = Tr(B_1 M) Q
  QLA_M_eq_c_times_M (&deriv, &trB1M, &Q);
  
  // deriv += Tr(B_2 M) Q^2
  QLA_M_peq_c_times_M (&deriv, &trB2M, &Q2);
  
  // deriv += f1 M
  QLA_M_peq_c_times_M (&deriv, &qf1, &matI);
  
  // deriv += f2 Q M
  QLA_M_eq_M_times_M  (&prod, &Q, &matI); // Q M
  QLA_M_peq_c_times_M (&deriv, &qf2, &prod); 
  
  // deriv += f2 M Q
  QLA_M_eq_M_times_M  (&prod, &matI, &Q); // M Q
  QLA_M_peq_c_times_M (&deriv, &qf2, &prod);

  QLA_M_eq_c_times_M  (&deriv, &minus_i, &deriv);
  
  printf("M = \n");
  printm(&matI);

  printf("deriv = \n");
  printm(&deriv);
  

  QLA_M_meq_M(&deriv, &expiQ);
  printf("diff = \n");
  printm(&deriv);
  /*
  printf("2/3f0 = \n");
  f0 *= 2.0/3;
  printc99(&f0);

  printc(&qf0);
  printc(&qf1);
  printc(&qf2);
  */

#endif

#if QDP_Colors == 2
  QLA_ColorMatrix expO;
  QLA_Complex Tr, det;
  QLA_c_eq_c_times_c(det, QLA_elem_M(O,0,0),QLA_elem_M(O,1,1));
  QLA_c_meq_c_times_c(det, QLA_elem_M(O,0,1), QLA_elem_M(O,1,0));

  QLA_C_eq_trace_M(&Tr, &O);
  QLA_Complex qs, qt;
  QLA_c_eq_r_times_c(qs, 0.5, Tr); // s=TrA/2

  QLA_Complex qs2;
  QLA_c_eq_c_times_c(qs2, qs, qs); //s2 = s^2
  QLA_c_meq_c(qs2, det); //s2 = s^2 - detA

  double _Complex t = QLA_real(qs2) + QLA_imag(qs2) * _Complex_I;
  t = csqrt(t); // sqrt(s^2 - det A)
  QLA_c_eq_r_plus_ir(qt, creal(t), cimag(t)); // t = sqrt(s^2 - det A)

  printf(" Matrix O = \n"); printm(&O);
  printf("TrO = ");         printc(&Tr);
  printf("detO = ");        printc(&det); 
  printf("s = ");           printc(&qs);
  printf("t^2 = ");         printc(&qs2);
  printf("t = ");           printc(&qt);

  //use the QLA exp function
  QLA_ColorMatrix qla_exp;
  QLA_M_eq_exp_M(&qla_exp, &O); //exp(O)
  
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

  //derivative of the exponential
  double _Complex f0s, f1s, f1t, f0t2, f1t2, t2;
  t2 = t*t;

  f0s = f0 - f1;
  f1s = f1;
  if (cabs(t) > 0.05) { 
    f1t = ((f0-f1) + s*f1)/t;
    f1t2 = f1t/t;
  }
  else { //when |t| < 0.05, the error is O(10^{-14})
    f1t = exps * t/3 * (1+t2/10*(1+t2/28));
    f1t2 = exps / 3 * (1+t2/10*(1+t2/28));
  }
 
  //  f0t = t * f1 - s * f1t;
  f0t2 = f1 - s * f1t2;
  
  printf("f0   = \n"); printc99(&f0);
  printf("f1   = \n"); printc99(&f1);
  printf("f0s  = \n"); printc99(&f0s);
  printf("f1s  = \n"); printc99(&f1s);

  printf("f1t  = \n"); printc99(&f1t);
  printf("f0t2 = \n"); printc99(&f0t2);
  printf("f1t2 = \n"); printc99(&f1t2);



  QLA_Complex qf0, qf1;
  QLA_c_eq_r_plus_ir(qf0, creal(f0), cimag(f0));
  QLA_c_eq_r_plus_ir(qf1, creal(f1), cimag(f1));

  QLA_M_eq_c_times_M(&expO, &qf1, &O);
  QLA_M_peq_c(&expO, &qf0);
  printf("QLA exp = \n"); printm(&qla_exp);
  printf("my expO = \n"); printm(&expO);

  /*
    QLA_Complex qf0s, qf0t, qf1s, qf1t;
    QLA_c_eq_r_plus_ir(qf0s, creal(f0s), cimag(f0s));
    QLA_c_eq_r_plus_ir(qf0t, creal(f0t), cimag(f0t));
    QLA_c_eq_r_plus_ir(qf1s, creal(f1s), cimag(f1s));
    QLA_c_eq_r_plus_ir(qf1t, creal(f1t), cimag(f1t));
  */

  QLA_ColorMatrix deriv;
  QLA_M_eq_zero(&deriv);

  QLA_ColorMatrix B, AB;
  QLA_M_eq_M(&B, &matI);

  //QLA_c_eq_r_plus_ir(QLA_elem_M(B,1,0), 0.1, 0.2);
  //QLA_c_eq_r_plus_ir(QLA_elem_M(B,0,1), 0.2, 0.1);
  
  printf("B=\n");  printm(&B);

  QLA_M_eq_M_times_M(&AB, &O, &B);
  printf("AB=\n"); printm(&AB);
  
  QLA_M_eq_c_times_M(&deriv, &qf1, &B); //f1 * B
  printf("f1 B = \n"); printm(&deriv);

  QLA_Complex trB, trAB;
  QLA_C_eq_trace_M(&trB,  &B); 
  QLA_C_eq_trace_M(&trAB, &AB);
  
  double _Complex ctrB  = QLA_real(trB) + _Complex_I * QLA_imag(trB);
  double _Complex ctrAB = QLA_real(trAB) + _Complex_I * QLA_imag(trAB);
  double _Complex coeff;
  
  coeff  = (f0s - f0t2 * s) * ctrB;
  coeff += (f1s - f1t2 * s) * ctrAB;
  coeff *= 0.5;
  
  printf("coeff = "); printc99(&coeff);
  QLA_Complex qc;
  QLA_D_c_eq_c99(qc, coeff);

  printc(&qc);

  QLA_M_peq_c_times_M(&deriv, &qc, &matI); // f1 * B + () I
  printf("f1B+()I=\n"); printm(&deriv);
  
  coeff = 0.5 * (f0t2 * ctrB + f1t2 * ctrAB);
  QLA_D_c_eq_c99(qc, coeff);

  printc(&qc);

  QLA_M_peq_c_times_M(&deriv, &qc, &O);

  printm(&deriv);

  exp_deriv_site(&deriv, &expO, &O, &B);
  printm(&deriv);

#endif

  return 0;
}
