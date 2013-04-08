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
}

static void
getfs(double _Complex *f0, double _Complex *f1, double _Complex *f2,
      double c0, double c1)
{
  double _Complex h0, h1, h2, e2iu, emiu;
  double c0m, t, u, w, u2, w2, cw, xi0, di;
  int sign=0;
  if(c0<0) { sign=1; c0=-c0; }
  c0m = 2.*pow(c1/3.,1.5);
  t = acos(c0/c0m);
  u = sqrt(c1/3.)*cos(t/3.);
  w = sqrt(c1)*sin(t/3.);
  u2 = u*u;
  w2 = w*w;
  cw = cos(w);
  if(w2>0.0025) xi0 = sin(w)/w;
  else xi0 = 1 - w2/6.*(1 - w2/20.*(1 - w2/42.));
  e2iu = cexp(2*I*u);
  emiu = cexp(-I*u);
  h0 = (u2-w2)*e2iu + emiu*(8*u2*cw+2*I*u*(3*u2+w2)*xi0);
  h1 = 2*u*e2iu - emiu*(2*u*cw-I*(3*u2-w2)*xi0);
  h2 = e2iu - emiu*(cw+3*I*u*xi0);
  di = 1/(9*u2-w2);
  *f0 = di*h0;
  *f1 = di*h1;
  *f2 = di*h2;
  if(sign) {
    *f0 = conj(*f0);
    *f1 = -conj(*f1);
    *f2 = conj(*f2);
  }
}

int
main(void)
{
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
  QLA_M_eq_Ma(&tmp, &O); // tmp = O^+;
  QLA_M_meq_M(&tmp, &O); // tmp = O^+ - O;
  QLA_M_eq_r_times_M(&iQ, &half, &tmp); // iQ = -(O^+ - O)/2;
  
  QLA_D_Real r = 0.5/QLA_Nc;
  QLA_M_eq_r_times_M(&tmp, &r, &tmp);
  QLA_C_eq_trace_M(&tr, &tmp); // trace iQ/2N

  printf("%f\n", r);
  printf("%f %f i\n", QLA_real(tr), QLA_imag(tr));
  
  QLA_M_meq_c_times_M(&iQ, &tr, &matI); //final iQ

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
