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
  QLA_ColorMatrix x, y, c, f, f2, yt;

  for(int i=0; i<QLA_Nc; i++) {
    for(int j=0; j<QLA_Nc; j++) {
      QLA_c_eq_r_plus_ir(QLA_elem_M(x,i,j), i+1, QLA_Nc*(j+1));
      QLA_c_eq_r_plus_ir(QLA_elem_M(c,i,j), QLA_Nc*(j+1), i-1);
      //QLA_c_eq_r(QLA_elem_M(x,i,j), 0);
      //QLA_c_eq_r(QLA_elem_M(c,i,j), 0);
    }
    QLA_c_eq_r_plus_ir(QLA_elem_M(x,i,i), 2+1, 1);
    QLA_c_eq_r_plus_ir(QLA_elem_M(c,i,i), 0.1, 0.5);
  }
  QLA_M_eq_sqrt_M(&y, &x);
  QLA_M_eq_Ma(&yt, &y);
  //sqrt_deriv_site(&f, &y, &x, &c);
  sylsolve_site(&f, &y, &y, &c);
  printm(&x);
  printm(&y);
  printm(&c);
  printm(&f);

  QLA_ColorMatrix dx, dy;
  for(int i=0; i<QLA_Nc; i++) {
    for(int j=0; j<QLA_Nc; j++) {
      QLA_Real eps = 1e-6;
      QLA_M_eq_M(&dx, &x);
      QLA_c_peq_r(QLA_elem_M(dx,i,j), eps);
      QLA_M_eq_sqrt_M(&dy, &dx);
      QLA_M_meq_M(&dy, &y);
      QLA_Real ieps = 1/eps;
      QLA_M_eq_r_times_M(&dy, &ieps, &dy);
      //printm(&dy);
      QLA_M_eq_M_times_Ma(&f2, &dy, &c);
      QLA_Real r;
      QLA_R_eq_re_trace_M(&r, &f2);
      printf("%g\n", r);
    }
  }
//  printm(&f2);
  return 0;
}
