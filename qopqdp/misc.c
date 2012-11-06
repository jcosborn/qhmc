#include "qhmc_qopqdp_common.h"
#include <string.h>
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
  double max = nrm;
  QMP_max_double(&max);
  return sqrt(max);
}

static int isfirst=1;
static void
setfirst(QLA_ColorMatrix *m, int i, void *args)
{
  QLA_ColorMatrix *r = (QLA_ColorMatrix*)args;
  if(isfirst) {
    isfirst = 0;
    memcpy(r, m, sizeof(QLA_ColorMatrix));
  }
}

static int badIndexPlus1=0;
static void
compare(QLA_ColorMatrix *m, int i, void *args)
{
  QLA_ColorMatrix *r = (QLA_ColorMatrix*)args;
  if(badIndexPlus1==0) {
    if(memcmp(m, r, sizeof(QLA_ColorMatrix))) badIndexPlus1 = i+1;
  }
}

int
check_uniform_M(QDP_ColorMatrix *m, QDP_Subset s)
{
  QLA_ColorMatrix r;
  isfirst = 1;
  QDP_M_eq_funcia(m, setfirst, &r, s);
  QMP_broadcast(&r, sizeof(QLA_ColorMatrix));
  badIndexPlus1 = 0;
  QDP_M_eq_funcia(m, compare, &r, s);
  return badIndexPlus1;
}
