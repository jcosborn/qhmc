#include "qhmc_qopqdp_common.h"
#include <math.h>

static void
infnorm(QLA_ColorMatrix *m, int i, void *args)
{
  QLA_Real *p = (QLA_Real*)args;
  QLA_Real r;
  QLA_r_eq_norm2_M(&r, m);
  if(r>*p) *p = r;
}

QLA_Real
infnorm_M(QDP_ColorMatrix *m, QDP_Subset s)
{
  QLA_Real nrm = 0;
  QDP_M_eq_funcia(m, infnorm, &nrm, s);
  return sqrt(nrm);
}
