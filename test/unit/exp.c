#include <qhmc_qopqdp_common.h>

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
}



int
main(void)
{
  QLA_ColorMatrix O, iQ, tmp, I;
  QLA_M_eq_zero(&I);

  for(int i=0; i<QLA_Nc; i++) {
    QLA_c_eq_r_plus_ir(QLA_elem_M(I,i,i), 1.0, 0.0);
  }
  printm(&I);

  QLA_Complex tr;
  QLA_Real half = 0.5;

  for(int i=0; i<QLA_Nc; i++) {
    for(int j=0; j<QLA_Nc; j++) {
      QLA_c_eq_r_plus_ir(QLA_elem_M(O,i,j), i+1, QLA_Nc*(j+1));
    }
    QLA_c_eq_r_plus_ir(QLA_elem_M(O,i,i), 2+1, 1);
  }
  QLA_M_eq_Ma(&tmp, &O); // tmp = O^+;
  QLA_M_meq_M(&tmp, &O); // tmp = O^+ - O;
  QLA_M_eq_r_times_M(&iQ, &half, &tmp); // iQ = -(O^+ - O)/2;
  
  QLA_D_Real r = 0.5/QLA_Nc;
  QLA_M_eq_r_times_M(&tmp, &r, &tmp);
  QLA_C_eq_trace_M(&tr, &tmp); // trace iQ/2N

  printf("%f\n", r);
  printf("%f %f i\n", QLA_real(tr), QLA_imag(tr));
  
  QLA_M_meq_c_times_M(&iQ, &tr, &I); //final iQ

  //use the QLA exp function
  QLA_ColorMatrix qla_exp;
  QLA_M_eq_exp_M(&qla_exp, &iQ); //exp(iQ)

  //use my own implementation
  QLA_ColorMatrix expiQ;
  
  

  printm(&O);
  printf("iQ = \n");
  printm(&iQ);
  printf("QLA: exp(iQ) = \n");
  printm(&qla_exp);
  
  return 0;
}
