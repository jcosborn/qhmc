#include <math.h>
#include "qhmc_qopqdp_common.h"
#ifdef HAVE_NC1
#include <qdp_f1.h>
#include <qop_f1.h>
#endif
#ifdef HAVE_NC2
#include <qdp_f2.h>
#include <qop_f2.h>
#endif
#ifdef HAVE_NC3
#include <qdp_f3.h>
#include <qop_f3.h>
#endif

static void
get_resid(QDP_ColorVector *r, QDP_ColorVector *x, QDP_ColorVector *b,
	  QOP_FermionLinksAsqtad *fla, QLA_Real m, QOP_evenodd_t eo, QDP_ColorVector *t)
{
  if(eo==QOP_EVEN) {
    QOP_asqtad_dslash_qdp(NULL, fla, m, r, x, QOP_ODD, QOP_EVEN);
    QOP_asqtad_dslash_qdp(NULL, fla, m, r, r, QOP_EVEN, QOP_ODD);
    QLA_Real mi = 1/m;
    QDP_V_eq_r_times_V_plus_V(r, &mi, r, b, QDP_even);
    QDP_V_meq_r_times_V(r, &m, x, QDP_even);
  } else if(eo==QOP_EVENODD) {
    QOP_asqtad_dslash_qdp(NULL, fla, m, r, x, QOP_EVENODD, QOP_EVENODD);
    QDP_V_eq_V_minus_V(r, b, r, QDP_all);
  } else {
    QOP_asqtad_dslash_qdp(NULL, fla, m, r, x, QOP_EVEN, QOP_ODD);
    QOP_asqtad_dslash_qdp(NULL, fla, m, r, r, QOP_ODD, QOP_EVEN);
    QLA_Real mi = 1/m;
    QDP_V_eq_r_times_V_plus_V(r, &mi, r, b, QDP_odd);
    QDP_V_meq_r_times_V(r, &m, x, QDP_odd);
  }
}

void
asqtadInvert(QOP_info_t *info, QOP_FermionLinksAsqtad *fla,
	     QOP_invert_arg_t *invarg, QOP_resid_arg_t *residarg[],
	     QLA_Real *masses, int nm, QDP_ColorVector *prop[],
	     QDP_ColorVector *source)
{
#define NC QDP_get_nc(source)
  QDP_Subset sub=QDP_all;
  if(invarg->evenodd==QOP_EVEN) {
    sub = QDP_even;
  } else if(invarg->evenodd==QOP_ODD) {
    sub = QDP_odd;
  }
  QOP_F_FermionLinksAsqtad *ffla = QOP_FD_asqtad_create_L_from_L(fla);
  QLA_F_Real m[nm];
  QDP_F_ColorVector *x[nm], *b;
  for(int i=0; i<nm; i++) {
    m[i] = masses[i];
    x[i] = QDP_F_create_V();
    QDP_F_V_eq_zero(x[i], sub);
    residarg[i]->final_iter = 0;
    residarg[i]->final_restart = 0;
  }
  b = QDP_F_create_V();
  QDP_FD_V_eq_V(b, source, sub);

  double flops=0, secs=0;
  int its=0;
  double preres = 1e-4;

  if(nm>1) {
    QOP_resid_arg_t **ra = residarg;
    QLA_F_Real *mp = m;
    QDP_F_ColorVector **xp = x;
    switch(QOP_Nc) {
#ifdef HAVE_NC1
    case 1: QOP_F1_asqtad_invert_multi_qdp(info, (QOP_F1_FermionLinksAsqtad*)ffla, invarg, &ra, &mp, &nm, (QDP_F1_ColorVector***)&xp, (QDP_F1_ColorVector**)&b, 1);
      break;
#endif
#ifdef HAVE_NC2
    case 2: QOP_F2_asqtad_invert_multi_qdp(info, (QOP_F2_FermionLinksAsqtad*)ffla, invarg, &ra, &mp, &nm, (QDP_F2_ColorVector***)&xp, (QDP_F2_ColorVector**)&b, 1);
      break;
#endif
#ifdef HAVE_NC3
    case 3: QOP_F3_asqtad_invert_multi_qdp(info, (QOP_F3_FermionLinksAsqtad*)ffla, invarg, &ra, &mp, &nm, (QDP_F3_ColorVector***)&xp, (QDP_F3_ColorVector**)&b, 1);
      break;
#endif
    default: QOP_F_asqtad_invert_multi_qdp(info, ffla, invarg, &ra, &mp, &nm, &xp, &b, 1);
    }
    for(int i=0; i<nm; i++) {
      QDP_DF_V_eq_V(prop[i], x[i], sub);
    }
    its += residarg[0]->final_iter;
    secs += info->final_sec;
    flops += info->final_flop;
  }

  QLA_Real norm2in;
  QDP_r_eq_norm2_V(&norm2in, source, sub);
  QDP_ColorVector *Dprop, *r;
  r = QDP_create_V();
  Dprop = QDP_create_V();

  for(int i=0; i<nm; i++) {
    QLA_Real rsq, rsqmin;
    int iters=0, restarts=0;
    QOP_resid_arg_t res_arg = *residarg[i];
    QLA_Real rsqstop = res_arg.rsqmin*norm2in;

    while(1) {
      get_resid(r, prop[i], source, fla, masses[i], invarg->evenodd, Dprop);
      QDP_r_eq_norm2_V(&rsq, r, sub);
      //printf0("its: %i  resid[%i]/in = %g\n", its, i, sqrt(rsq/norm2in));
      if(rsq<rsqstop) break;
      if(restarts>invarg->max_restarts) break;
      if(iters>0) restarts++;

      rsqmin = preres*preres;
      if(rsqstop/rsq>rsqmin) rsqmin = rsqstop/rsq;
      rsqmin *= 0.99;
      res_arg.rsqmin = rsqmin;

      if(1) {
	QDP_FD_V_eq_V(b, r, sub);
	QDP_F_V_eq_zero(x[i], sub);
	switch(QOP_Nc) {
#ifdef HAVE_NC1
	case 1: QOP_F1_asqtad_invert_qdp(info, (QOP_F1_FermionLinksAsqtad*)ffla, invarg, &res_arg, m[i], (QDP_F1_ColorVector*)x[i], (QDP_F1_ColorVector*)b);
	  break;
#endif
#ifdef HAVE_NC2
	case 2: QOP_F2_asqtad_invert_qdp(info, (QOP_F2_FermionLinksAsqtad*)ffla, invarg, &res_arg, m[i], (QDP_F2_ColorVector*)x[i], (QDP_F2_ColorVector*)b);
	  break;
#endif
#ifdef HAVE_NC3
	case 3: QOP_F3_asqtad_invert_qdp(info, (QOP_F3_FermionLinksAsqtad*)ffla, invarg, &res_arg, m[i], (QDP_F3_ColorVector*)x[i], (QDP_F3_ColorVector*)b);
	  break;
#endif
	default: QOP_F_asqtad_invert_qdp(info, ffla, invarg, &res_arg, m[i], x[i], b);
	}
	QDP_DF_V_eq_V(Dprop, x[i], sub);
      } else {
	QDP_V_eq_zero(Dprop, sub);
	QOP_asqtad_invert_qdp(info, fla, invarg, &res_arg, masses[i], Dprop, r);
      }
      residarg[i]->final_restart++;
      residarg[i]->final_iter += res_arg.final_iter;
      its += res_arg.final_iter;
      iters += res_arg.final_iter;
      secs += info->final_sec;
      flops += info->final_flop;
      QDP_V_peq_V(prop[i], Dprop, sub);
    }
    //printf0("iters = %4i  secs = %8f  mflops = %.0f resid = %g\n", its,
    //secs, 1e-6*flops/secs, sqrt(rsq/norm2in));
  }
  info->final_sec = secs;
  info->final_flop = flops;
  QDP_destroy_V(r);
  QDP_destroy_V(Dprop);
  QDP_F_destroy_V(b);
  for(int i=0; i<nm; i++) {
    QDP_F_destroy_V(x[i]);
  }
  QOP_F_asqtad_destroy_L(ffla);
#undef NC
}
