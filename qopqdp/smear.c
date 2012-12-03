#include "qhmc_qopqdp_common.h"

static void
zero_asqtad_coeffs(QOP_asqtad_coeffs_t *coeffs)
{
  coeffs->one_link = 0;
  coeffs->three_staple = 0;
  coeffs->five_staple = 0;
  coeffs->seven_staple = 0;
  coeffs->lepage = 0;
  coeffs->naik = 0;
}

static void
get_asqtad_coeffs(lua_State *L, int idx, QOP_asqtad_coeffs_t *coeffs)
{
  int cidx = lua_absindex(L, idx);
  tableLoopKeys(L, cidx) {
    //printf0("%s\n", tableLoopKey);
    tableLoopKeysIfKeySetDouble(L, "one_link", coeffs->one_link);
    else tableLoopKeysIfKeySetDouble(L, "three_staple", coeffs->three_staple);
    else tableLoopKeysIfKeySetDouble(L, "five_staple", coeffs->five_staple);
    else tableLoopKeysIfKeySetDouble(L, "seven_staple", coeffs->seven_staple);
    else tableLoopKeysIfKeySetDouble(L, "lepage", coeffs->lepage);
    tableLoopKeysEnd(L);
  }
}

static int
qopqdp_smear(lua_State *L)
{
  qassert(lua_gettop(L)==3);
  int nsg; get_table_len(L, 1, &nsg);
  gauge_t *sg[nsg]; qopqdp_gauge_array_check(L, 1, nsg, sg);
  int ng; get_table_len(L, 2, &ng);
  gauge_t *g[ng]; qopqdp_gauge_array_check(L, 2, ng, g);
  const char *type = tableGetString(L, 3, "type");
  if(strcmp(type,"fat7")==0) {
    qassert(nsg==1 && ng==1);
    QOP_asqtad_coeffs_t coeffs;
    zero_asqtad_coeffs(&coeffs);
    tableGetField(L, 3, "coeffs");
    get_asqtad_coeffs(L, -1, &coeffs);
    lua_pop(L, 1);
    coeffs.naik = 0;
    QOP_info_t info;
    //printf0("one_link = %g\n", coeffs.one_link);
    QOP_smear_fat7l_qdp(&info, sg[0]->links, g[0]->links, &coeffs);
    //for(int mu=0; mu<4; mu++) QDP_M_eq_M(sg[0]->links[mu], g[0]->links[mu], QDP_all);
  } else
  if(strcmp(type,"product")==0) {
    qassert(nsg==1 && ng==2);
    tableGetField(L, 3, "adj");
    int na; get_table_len(L, -1, &na);
    int adj[na]; get_bool_array(L, -1, na, adj);
    lua_pop(L, 1);
    qassert(na==ng);
    int nd = sg[0]->nd;
    if(adj[0]) {
      if(adj[1]) {
	for(int mu=0; mu<nd; mu++) {
	  QDP_M_eq_Ma_times_Ma(sg[0]->links[mu], g[0]->links[mu], g[1]->links[mu], QDP_all);
	}
      } else {
	for(int mu=0; mu<nd; mu++) {
	  QDP_M_eq_Ma_times_M(sg[0]->links[mu], g[0]->links[mu], g[1]->links[mu], QDP_all);
	}
      }
    } else {
      if(adj[1]) {
	for(int mu=0; mu<nd; mu++) {
	  QDP_M_eq_M_times_Ma(sg[0]->links[mu], g[0]->links[mu], g[1]->links[mu], QDP_all);
	}
      } else {
	for(int mu=0; mu<nd; mu++) {
	  QDP_M_eq_M_times_M(sg[0]->links[mu], g[0]->links[mu], g[1]->links[mu], QDP_all);
	}
      }
    }
  } else
  if(strcmp(type,"antiherm")==0) {
    qassert(nsg==1 && ng==1);
    int nd = sg[0]->nd;
    for(int mu=0; mu<nd; mu++) {
      //QDP_M_eq_M(sg[0]->links[mu], g[0]->links[mu], QDP_all);
      //QDP_M_eq_antiherm_M(sg[0]->links[mu], g[0]->links[mu], QDP_all);
      QDP_M_eq_M(sg[0]->links[mu], g[0]->links[mu], QDP_all);
      QDP_M_meq_Ma(sg[0]->links[mu], g[0]->links[mu], QDP_all);
      QLA_Real half=0.5;
      QDP_M_eq_r_times_M(sg[0]->links[mu], &half, sg[0]->links[mu], QDP_all);
    }
  } else
  if(strcmp(type,"tracelessAntiherm")==0) {
    qassert(nsg==1 && ng==1);
    int nd = sg[0]->nd;
    for(int mu=0; mu<nd; mu++) {
      //QDP_M_eq_M(sg[0]->links[mu], g[0]->links[mu], QDP_all);
      QDP_M_eq_antiherm_M(sg[0]->links[mu], g[0]->links[mu], QDP_all);
    }
  } else
  if(strcmp(type,"mobius")==0) {
    qassert(nsg==1 && ng==1);
    tableGetField(L, 3, "coeffs");
    int nc; get_table_len(L, -1, &nc);
    qassert(nc==4);
    double coeffs[nc]; get_double_array(L, -1, nc, coeffs);
    lua_pop(L, 1);
    QLA_Real a=coeffs[0], b=coeffs[1], c=coeffs[2], d=coeffs[3];
    int nd = sg[0]->nd;
    for(int mu=0; mu<nd; mu++) {
      int i;
      QDP_loop_sites(i, QDP_all, ({
	  QLA_ColorMatrix t1, t2, t3;
	  QLA_ColorMatrix *gmui = QDP_site_ptr_readonly_M(g[0]->links[mu], i);
	  QLA_ColorMatrix *sgmui= QDP_site_ptr_readwrite_M(sg[0]->links[mu],i);
	  QLA_M_eq_r_times_M(&t1, &b, gmui);
	  QLA_M_eq_r_times_M(&t2, &d, gmui);
	  for(int ic=0; ic<QLA_Nc; ic++) {
	    QLA_c_peq_r(QLA_elem_M(t1,ic,ic), a);
	    QLA_c_peq_r(QLA_elem_M(t2,ic,ic), c);
	  }
	  QLA_M_eq_inverse_M(&t3, &t2);
	  QLA_M_eq_M_times_M(sgmui, &t1, &t3);
	  }));
    }
  } else
  if(strcmp(type,"exp")==0) {
    qassert(nsg==1 && ng==1);
    tableGetField(L, 3, "rho");
    int ns; get_table_len(L, -1, &ns);
    double rho[ns]; get_double_array(L, -1, ns, rho);
    lua_pop(L, 1);
#if 1
    int nd = sg[0]->nd;
    for(int mu=0; mu<nd; mu++) {
      QLA_Real r = rho[mu];
      QDP_M_eq_r_times_M(sg[0]->links[mu], &r, g[0]->links[mu], QDP_all);
    }
#else    
    QDP_ColorMatrix *cm = QDP_create_M();
    for(int mu=0; mu<ns; mu++) {
      QLA_Real r = rho[mu];
      if(r) {
	QDP_M_eq_r_times_M(cm, &r, g[0]->links[mu], QDP_all);
	QDP_M_eq_exp_M(sg[0]->links[mu], cm, QDP_all);
      } else {
	QLA_Real one = 1;
	QDP_M_eq_r(sg[0]->links[mu], &one, QDP_all);
      }
    }
    QDP_destroy_M(cm);
#endif
  } else {
    printerr("unknown smearing type: %s\n", type);
    ABORT(1);
  }
  return 0;
}

static int
qopqdp_smearChain(lua_State *L)
{
  qassert(lua_gettop(L)==5);
  int nf; get_table_len(L, 1, &nf);
  force_t *f[nf]; qopqdp_force_array_check(L, 1, nf, f);
  int nfc; get_table_len(L, 2, &nfc);
  force_t *fc[nfc]; qopqdp_force_array_check(L, 2, nfc, fc);
  int nsg; get_table_len(L, 3, &nsg);
  gauge_t *sg[nsg]; qopqdp_gauge_array_check(L, 3, nsg, sg);
  int ng; get_table_len(L, 4, &ng);
  gauge_t *g[ng]; qopqdp_gauge_array_check(L, 4, ng, g);
  const char *type = tableGetString(L, 5, "type");
  if(strcmp(type,"fat7")==0) {
    qassert(nf==1 && nfc==1 && nsg==1 && ng==1);
    QOP_asqtad_coeffs_t coeffs;
    zero_asqtad_coeffs(&coeffs);
    tableGetField(L, 5, "coeffs");
    get_asqtad_coeffs(L, -1, &coeffs);
    lua_pop(L, 1);
    coeffs.naik = 0;
    //printf0("one_link = %g\n", coeffs.one_link);
    QOP_info_t info;
    //for(int mu=0; mu<4; mu++) QDP_M_eq_zero(f[0]->force[mu], QDP_all);
    //for(int mu=0; mu<4; mu++) QDP_M_peq_M(f[0]->force[mu], fc[0]->force[mu], QDP_all);
    //coeffs.one_link *= 0.5;
    //coeffs.three_staple *= -1;
    QOP_asqtad_deriv(&info, g[0]->links, f[0]->force, &coeffs,
		     fc[0]->force, NULL);
  } else
  if(strcmp(type,"product")==0) {
    qassert(nf==2 && nfc==1 && nsg==1 && ng==2);
    tableGetField(L, 5, "adj");
    int na; get_table_len(L, -1, &na);
    int adj[na]; get_bool_array(L, -1, na, adj);
    lua_pop(L, 1);
    qassert(na==ng);
    int nd = f[0]->nd;
    if(adj[0]) {
      if(adj[1]) {
	for(int mu=0; mu<nd; mu++) {
	  QDP_M_peq_Ma_times_Ma(f[0]->force[mu], g[1]->links[mu], fc[0]->force[mu], QDP_all);
	  QDP_M_peq_Ma_times_Ma(f[1]->force[mu], fc[0]->force[mu], g[0]->links[mu], QDP_all);
	}
      } else {
	for(int mu=0; mu<nd; mu++) {
	  QDP_M_peq_M_times_Ma(f[0]->force[mu], g[1]->links[mu], fc[0]->force[mu], QDP_all);
	  QDP_M_peq_M_times_M(f[1]->force[mu], g[0]->links[mu], fc[0]->force[mu], QDP_all);
	}
      }
    } else {
      if(adj[1]) {
	for(int mu=0; mu<nd; mu++) {
	  QDP_M_peq_M_times_M(f[0]->force[mu], fc[0]->force[mu], g[1]->links[mu], QDP_all);
	  QDP_M_peq_Ma_times_M(f[1]->force[mu], fc[0]->force[mu], g[0]->links[mu], QDP_all);
	}
      } else {
	for(int mu=0; mu<nd; mu++) {
	  QDP_M_peq_M_times_Ma(f[0]->force[mu], fc[0]->force[mu], g[1]->links[mu], QDP_all);
	  QDP_M_peq_Ma_times_M(f[1]->force[mu], g[0]->links[mu], fc[0]->force[mu], QDP_all);
	}
      }
    }
  } else
  if(strcmp(type,"antiherm")==0) {
    qassert(nf==1 && nfc==1 && nsg==1 && ng==1);
    int nd = f[0]->nd;
    for(int mu=0; mu<nd; mu++) {
      /// FIXME: make it add to force
      QDP_M_eq_M(f[0]->force[mu], fc[0]->force[mu], QDP_all);
      QDP_M_meq_Ma(f[0]->force[mu], fc[0]->force[mu], QDP_all);
      QLA_Real half=0.5;
      QDP_M_eq_r_times_M(f[0]->force[mu], &half, f[0]->force[mu], QDP_all);
    }
  } else
  if(strcmp(type,"tracelessAntiherm")==0) {
    // X = 0.5*(A-A^+)
    // Y = X - <X>/Nc
    // <Y C^+ + Y^+ C> = <X C^+ + X^+ C> - (<X><C^+>+<X^+><C>)/Nc
    // = <X(C^+ - C)> - <X>(<C^+>-<C>)/Nc
    // dA -> 0.5*(C-C^+) 
    qassert(nf==1 && nfc==1 && nsg==1 && ng==1);
    int nd = f[0]->nd;
    for(int mu=0; mu<nd; mu++) {
      //QDP_M_eq_M(f[0]->force[mu], fc[0]->force[mu], QDP_all);
      /// FIXME: make it add to force
      QDP_M_eq_antiherm_M(f[0]->force[mu], fc[0]->force[mu], QDP_all);
    }
  } else
  if(strcmp(type,"mobius")==0) {
    qassert(nf==1 && nfc==1 && nsg==1 && ng==1);
    tableGetField(L, 5, "coeffs");
    int nc; get_table_len(L, -1, &nc);
    qassert(nc==4);
    double coeffs[nc]; get_double_array(L, -1, nc, coeffs);
    lua_pop(L, 1);
    QLA_Real a=coeffs[0], b=coeffs[1], c=coeffs[2], d=coeffs[3];
    QLA_Real detinv = -1/(a*d-b*c);
    int nd = f[0]->nd;
    for(int mu=0; mu<nd; mu++) {
      int i;
      QDP_loop_sites(i, QDP_all, ({
	  QLA_ColorMatrix t1, t2, t3;
	  QLA_ColorMatrix *fmui = QDP_site_ptr_readwrite_M(f[0]->force[mu], i);
	  QLA_ColorMatrix *fcmui = QDP_site_ptr_readonly_M(fc[0]->force[mu],i);
	  QLA_ColorMatrix *sgmui = QDP_site_ptr_readonly_M(sg[0]->links[mu],i);
	  QLA_M_eq_r_times_M(&t1, &d, sgmui);
	  for(int ic=0; ic<QLA_Nc; ic++) {
	    QLA_c_meq_r(QLA_elem_M(t1,ic,ic), b);
	  }
	  QLA_M_eq_Ma_times_M(&t2, &t1, fcmui);
	  QLA_M_eq_M_times_Ma(&t3, &t2, &t1);
	  QLA_M_peq_r_times_M(fmui, &detinv, &t3);
	  }));
    }
  } else
  if(strcmp(type,"exp")==0) {
    qassert(nf==1 && nfc==1 && nsg==1 && ng==1);
    tableGetField(L, 5, "rho");
    int ns; get_table_len(L, -1, &ns);
    double rho[ns]; get_double_array(L, -1, ns, rho);
    lua_pop(L, 1);
#if 1
    int nd = sg[0]->nd;
    for(int mu=0; mu<nd; mu++) {
      QLA_Real r = rho[mu];
      QDP_M_peq_r_times_M(f[0]->force[mu], &r, fc[0]->force[mu], QDP_all);
    }
#else    
    QDP_ColorMatrix *cm2 = QDP_create_M();
    for(int mu=0; mu<ns; mu++) {
      QLA_Real r = rho[mu];
      if(r) {
	QDP_ColorMatrix *cm1 = sg[0]->links[mu];
	QDP_M_eq_M_times_Ma(cm1, g[1]->links[mu], g[0]->links[mu], QDP_all);
	QDP_M_eq_antiherm_M(cm2, cm1, QDP_all);
	QDP_M_eq_r_times_M(cm1, &r, cm2, QDP_all);
	QDP_M_eq_exp_M(cm2, cm1, QDP_all);
	QDP_M_eq_M_times_M(sg[0]->links[mu], cm2, g[0]->links[mu], QDP_all);
      }
    }
    QDP_destroy_M(cm2);
#endif
  } else {
    printerr("unknown smearing type: %s\n", type);
    ABORT(1);
  }
  return 0;
}

static struct luaL_Reg qopqdp_reg[] = {
  { "smear",      qopqdp_smear },
  { "smearChain", qopqdp_smearChain },
  { NULL, NULL}
};

void
open_qopqdp_smear(lua_State* L)
{
  luaL_register(L, "qopqdp", qopqdp_reg);
}
