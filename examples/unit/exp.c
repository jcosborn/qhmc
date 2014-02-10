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

int
main(void) {
  QLA_ColorMatrix O, iQ, tmp, matI;
  QLA_M_eq_zero(&matI);

  for(int i=0; i<QLA_Nc; i++) {
    QLA_c_eq_r_plus_ir(QLA_elem_M(matI,i,i), 1.0, 0.0);
  }
  printm(&matI);

  QLA_Complex tr;
  QLA_Real half = 0.5;

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
  QLA_Complex s, t;
  QLA_c_eq_r_times_c(s, 0.5, Tr); // s=TrA/2

  QLA_Complex s2;
  QLA_c_eq_c_times_c(s2, s, s); //s2 = s^2
  QLA_c_meq_c(s2, det); //s2 = s^2 - detA

  double _Complex dc_t = QLA_real(s2) + QLA_imag(s2) * _Complex_I;
  dc_t = csqrt(dc_t); // sqrt(s^2 - det A)
  QLA_c_eq_r_plus_ir(t, creal(dc_t), cimag(dc_t)); // t = sqrt(s^2 - det A)

  printf(" Matrix O = \n"); printm(&O);
  printf("TrO = ");         printc(&Tr);
  printf("detO = ");        printc(&det); 
  printf("s = ");           printc(&s);
  printf("t^2 = ");         printc(&s2);
  printf("t = ");           printc(&t);

  //use the QLA exp function
  QLA_ColorMatrix qla_exp;
  QLA_M_eq_exp_M(&qla_exp, &O); //exp(O)
  
  double _Complex cosht, sinht, sinht_t;

  if(QLA_real(t) == 0 && QLA_imag(t) == 0) {
    cosht = 1;
    sinht = 0;
    sinht_t = 1;	
  } else {
    cosht = ccosh(dc_t);
    sinht = csinh(dc_t);
    sinht_t = sinht/dc_t;
  }

  double _Complex dc_s = QLA_real(s) + QLA_imag(s) * _Complex_I;
  double _Complex dc_f0, dc_f1;
  
  dc_f0 = cexp(dc_s) * (cosht - dc_s * sinht_t);
  dc_f1 = cexp(dc_s) * sinht_t;

  QLA_Complex f0, f1;
  QLA_c_eq_r_plus_ir(f0, creal(dc_f0), cimag(dc_f0));
  QLA_c_eq_r_plus_ir(f1, creal(dc_f1), cimag(dc_f1));

  QLA_M_eq_c_times_M(&expO, &f1, &O);
  QLA_M_peq_c(&expO, &f0);
  printf("QLA exp = \n"); printm(&qla_exp);
  printf("my expO = \n"); printm(&expO);


#endif

  return 0;
}
