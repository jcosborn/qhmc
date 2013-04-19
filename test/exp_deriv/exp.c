#include <qhmc_qopqdp_common.h>
#include <complex.h>
#include <math.h>

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

#if QDP_Colors == 3
  QLA_Complex ci;
  QLA_c_eq_r_plus_ir(ci, 0, 1);

  //use my own implementation
  QLA_ColorMatrix expiQ;
  QLA_ColorMatrix QQ;
  QLA_ColorMatrix QQQ;
  QLA_M_eq_M_times_M(&QQ, &iQ, &iQ); //-Q^2
  QLA_M_eq_M_times_M(&QQQ, &QQ, &iQ); //-iQ^3
  QLA_M_eq_c_times_M(&QQQ, &ci, &QQQ); //Q^3
 
  QLA_Complex c0, c1;
  QLA_C_eq_trace_M(&c0, &QQQ);
  QLA_c_eq_r_times_c(c0, 0.3333333, c0);

  QLA_C_eq_trace_M(&c1, &QQ);
  QLA_c_eq_r_times_c(c1, -0.5, c1);

  double _Complex tf0, tf1, tf2;
  getfs(&tf0, &tf1, &tf2, QLA_real(c0), QLA_real(c1));

  QLA_Complex f0, f1, f2;
  f0 = tf0; 
  f1 = tf1;
  f2 = tf2;

  printm(&O);
  printf("iQ = \n");
  printm(&iQ);
  printf("QLA: exp(iQ) = \n");
  printm(&qla_exp);
  printf("Q^3 = \n");
  printm(&QQQ);
  printf("c0 = "); printc(&c0);
  printf("c1 = "); printc(&c1);  
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

  
#endif

  return 0;
}
