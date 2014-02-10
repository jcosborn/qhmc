#include "qhmc_qopqdp_common.h"
#include <string.h>
#include <math.h>

// This is a private function to extract the SU(2) submatrix V from an SU(3) matrix.
// Taken from Chroma function "su2Extract_t". The submatrix correspond with sunFill.
#define NC nc
static void
su2Extract_local(NCPROT QLA_Real* r, QLA_ColorMatrix* v, int su2_index)
{
	// Determine the SU(3) indices corresponding to the SU(2) indices
	// of the SU(2) subgroup 'su2_index'
	int i1=0, i2;
	int found = 0;
	int del_i = 0;
	int index = -1;
	
	while (del_i < (QLA_Nc - 1) && found == 0)
	{
		del_i++;
		for (i1=0; i1 < (QLA_Nc - del_i); i1++)
		{
			index++;
			if (index == su2_index)
			{
				found = 1;
				break;
			}
		}
	}
	i2 = i1 + del_i;
	
	qassert(found != 0);

	// Compute the b(k) of A_SU(2) = b0 + i sum_k bk sigma_k	
	QLA_Complex s11, s12, s21, s22;
	//QLA_Real thing1, thing2;
	
	QLA_C_eq_elem_M(&s11, v, i1, i1);
	QLA_C_eq_elem_M(&s12, v, i1, i2);
	QLA_C_eq_elem_M(&s21, v, i2, i1);
	QLA_C_eq_elem_M(&s22, v, i2, i2);
	
	r[0] = QLA_real(s11) + QLA_real(s22);
	r[1] = QLA_imag(s12) + QLA_imag(s21);
	r[2] = QLA_real(s12) - QLA_real(s21);
	r[3] = QLA_imag(s11) - QLA_imag(s22);
	
	return;
}
#undef NC

// This is a private function to fill an SU(3) matrix V with a submatrix.
// Taken from Chroma function "sunFill", specified for SU(3).
#define NC nc
static void
sunFill_local(NCPROT QLA_ColorMatrix* v, QLA_Real* a, int su2_index)
{
	// Determine the SU(3) indices corresponding to the SU(2) indices
	// of the SU(2) subgroup 'su2_index'
	int i1, i2;
	int found = 0;
	int del_i = 0;
	int index = -1;
	
	while (del_i < (QLA_Nc - 1) && found == 0)
	{
		del_i++;
		for (i1 = 0; i1 < (QLA_Nc - del_i); i1++)
		{
			index++;
			if (index == su2_index)
			{
				found = 1;
				break;
			}
		}
	}
	i2 = i1 + del_i;
	
	qassert(found != 0);
	
	// Insert the b[k] of A_SU(2) = b0 + i sum_k bk sigma_k
	// back into the SU(N) matrix.
	
	
	QLA_Complex one1;
	QLA_Real re = 1.0;
	QLA_Real im = 0.0;
	QLA_ColorMatrix(v2);
	QLA_M_eq_zero(&v2);
	QLA_C_eq_R_plus_i_R(&one1, &re, &im);
	QLA_M_eq_elem_C(&v2, &one1, 1, 1);
	QLA_M_eq_elem_C(&v2, &one1, 2, 2);
	QLA_M_eq_elem_C(&v2, &one1, 0, 0);
	
	QLA_Complex fills;
	
	QLA_real(fills) = a[0];
	QLA_imag(fills) = a[3];
	QLA_M_eq_elem_C(&v2, &fills, i1, i1);
	
	QLA_real(fills) = a[2];
	QLA_imag(fills) = a[1];
	QLA_M_eq_elem_C(&v2, &fills, i1, i2);
	
	QLA_real(fills) = -a[2];
	QLA_imag(fills) = a[1];
	QLA_M_eq_elem_C(&v2, &fills, i2, i1);
	
	QLA_real(fills) = a[0];
	QLA_imag(fills) = -a[3];
	QLA_M_eq_elem_C(&v2, &fills, i2, i2);
	
	QLA_M_eq_M(v, &v2);
	
	return;
}
#undef NC

struct grelax_params
{
	int su2_index;
	int OrDo;
	double OrPara;
	int nd;
};

static double fuzz = 1e-12;

// Local function to perform the relaxation on the SU(2) subcomponents. 
#define NC nc
static void
private_local_grelax(NCPROT QLA_ColorMatrix(*v), int site, void *args)
{
	struct grelax_params* params = (struct grelax_params*)(args);
	QLA_Real r[4];
	QLA_Real a[4];
	QLA_Real r_l;
	QLA_Real lftmp;
	QLA_ColorMatrix(v2);
	
	QLA_M_eq_M(&v2, v);
	su2Extract_local(NCARG r, &v2, params->su2_index);
	
	// Now project onto SU(2)
	// r_l = sqrt(r[0]*r[0] + ... + r[3]*r[3])
	
	r_l = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);
	
	// Normalize, with a check to make sure things aren't too small.
	// Fill (r[0]/r_l, -r[1]/r_l, -r[2]/r_l, -r[3]/r_l) for r_l > fuzz
	// and (1,0,0,0) elsewhere.
	if (r_l > fuzz)
	{
		lftmp = 1.0/r_l;
		a[0] = r[0]*lftmp;
		a[1] = -r[1]*lftmp;
		a[2] = -r[2]*lftmp;
		a[3] = -r[3]*lftmp;
	}
	else
	{
		a[0] = 1.0;
		a[1] = a[2] = a[3] = 0.0;
	}
	

	
	// Now do the overrelaxation, if desired.
	if (params->OrDo > 0)
	{
		QLA_Real theta_old, oldsin, theta_new;
		QLA_Real overrelax_param = params->OrPara;
		
		// Get angle.
		theta_old = acos(a[0]);
		
		// Old sin.
		oldsin = sin(theta_old);
		
		// Overrelax, i.e., multiply by the angle.
		theta_new = overrelax_param*theta_old;
		
		// Compute sin(new)/sin(old).
		// Set the ratio to 0 if sin(old) < FUZZ.
		if (oldsin < fuzz)
		{
			lftmp = 0;
		}
		else
		{
			lftmp = sin(theta_new) / oldsin;
		}
		
		// get the new cos = a[0]
		a[0] = cos(theta_new);
		
		// get the new a_k, k = 1, 2, 3
		a[1] *= lftmp;
		a[2] *= lftmp;
		a[3] *= lftmp;
	}
	
	// Now fill the SU(3) matrix V with the SU(2) submatrix 'su2_index'
	// parameterized by a_k in the sigma matrix basis.
	
	sunFill_local(NCARG &v2, a, params->su2_index);
	QLA_M_eq_M(v, &v2);

}
#undef NC

// This is a private gauge relaxation function for gauge fixing.
// Based directly on the Chroma function of the name "grelax"
static void
private_grelax(QDP_ColorMatrix* g, gauge_t* u, int j_decay, int su2_index,
	       int cb, int OrDo, double OrPara)
{
#define NC QDP_get_nc(g)
	QDP_ColorMatrix *v, *extra, *extra2;
	QDP_ColorMatrix **ug, **u_neg;
	QDP_Real** r;
	QDP_Real** a;
	QDP_Real *r_l, *r_l2;
	//QDP_Real *theta_old, *oldsin, *theta_new;
	//QLA_Real overrelax_param;
	int nd = QDP_ndim();
	int mu;
	QDP_Subset checkboard;
	struct grelax_params subgroup_params;

	r = malloc(nd*sizeof(QDP_Real*));
	a = malloc(nd*sizeof(QDP_Real*));
	ug = malloc(nd*sizeof(QDP_ColorMatrix*));
	u_neg = malloc(nd*sizeof(QDP_ColorMatrix*));
	
	v = QDP_create_M();
	extra = QDP_create_M();
	extra2 = QDP_create_M();
	r_l = QDP_create_R();
	r_l2 = QDP_create_R();
	for (mu = 0; mu < nd; mu++)
	{
		r[mu] = QDP_create_R();
		a[mu] = QDP_create_R();
		ug[mu] = QDP_create_M();
		u_neg[mu] = QDP_create_M();
	}
	
	// Specify a checkerboard.
	if (cb == 0) // Even
	{
		checkboard = QDP_even;
	}
	else // Odd
	{
		checkboard = QDP_odd;
	}
	
	// Other things.
	int tDir = nd - 1;

	// ESW note: We may need to interchange the forward and back...
	// Rotate the matrices using the current gauge rotation.
	for (mu = 0; mu < nd; mu++)
	{
		// First get g(x+mu).
		QDP_M_eq_sM(ug[mu], g, QDP_neighbor[mu], QDP_forward, QDP_all);
		// Next get u[mu](x)*g(x+mu)^dagger.
		QDP_M_eq_M_times_Ma(extra, u->links[mu], ug[mu], QDP_all);
		// Lastly, get g(x)*u[mu](x)*g(x+mu)^dagger.
		QDP_M_eq_M_times_M(ug[mu], g, extra, QDP_all);
	}
	
	// Gather the Nd negative links attached to a site:
	// u_neg(x,mu) = U(x-mu, mu)
	for (mu = 0; mu < nd; mu++)
	{
		QDP_M_eq_sM(u_neg[mu], ug[mu], QDP_neighbor[mu], QDP_backward, QDP_all);
	}
	
	// Sum links to be gauge transformed on site x not in the direction
	// j_decay into matrix V:
	
	QDP_M_eq_zero(v, QDP_all);
	// This is split off because of Chroma's aniso support.
	if (tDir != j_decay)
	{
		// v[cb] = ug[tdir] + adj(u_neg[tdir])
		QDP_M_eq_Ma(extra, u_neg[tDir], checkboard);
		QDP_M_eq_M_plus_M(v, ug[tDir], extra, checkboard);
	}
	
	for (mu = 0; mu < nd; ++mu)
	{
		if (mu != j_decay && mu != tDir)
		{
			// v[cb] += ug[mu] + adj(u_neg[mu])
			QDP_M_eq_Ma(extra, u_neg[mu], checkboard);
			QDP_M_eq_M_plus_M(extra, ug[mu], extra, checkboard);
			QDP_M_eq_M(extra2, v, checkboard);
			QDP_M_eq_M_plus_M(v, extra2, extra, checkboard);
		}
	}
	
	// Once we have V, we can do everything to V on a local level.
	subgroup_params.su2_index = su2_index;
	subgroup_params.nd = nd;
	subgroup_params.OrDo = OrDo;
	subgroup_params.OrPara = OrPara;
	
	QDP_M_eq_funcia(v, private_local_grelax, &subgroup_params, checkboard);
	
	// End local v level stuff.
	
	// Now do the gauge transformation with matrix V:
	// Multiply into the global 'g' field
	QDP_M_eq_M_times_M(extra, v, g, checkboard);
	QDP_M_eq_M(g, extra, checkboard);
	
	// Take care of things, clean up memory, etc. 
	QDP_destroy_M(v);
	QDP_destroy_M(extra);
	QDP_destroy_M(extra2);
	QDP_destroy_R(r_l);
	QDP_destroy_R(r_l2);
	//QDP_destroy_I(lbtmp);
	//QDP_destroy_R(fuzz);
	//QDP_destroy_R(lftmp);
	for (mu = 0; mu < nd; mu++)
	{
		QDP_destroy_R(r[mu]);
		QDP_destroy_R(a[mu]);
		QDP_destroy_M(ug[mu]);
		QDP_destroy_M(u_neg[mu]);
	}
	
	free(r); free(a); free(ug); free(u_neg);
	
	return;
#undef NC
}

// This function performs gauge fixing to Coulomb or Landau gauge.
// It's based directly on Chroma's code, defined in coulgauge.cc.
// Arguments:
// 1. Gauge field, overrides it. Make a copy in advance if you want it.
// 2. Direction to fix orthogonally to. If this is > number dimensions, it'll fix to Landau.
//    Assume that x = 0, y = 1, z = 2, t = 3.
// 3. GF Accuracy -- accuracy of gauge fixing.
// 4. GF Max -- maximum number of iterations.
// 5. Use overrelaxation? (1 for yes, 0 for no) [OPTIONAL!]
// 6. Overrelaxation parameter [OPTIONAL!]
// Note that if you specify 5, you also have to specify 6. Default behavior
// is to not use overrelaxation.
// THIS CODE ASSUMES AN ISOTROPIC, SU(3) LATTICE.
int
qopqdp_gauge_coulomb(lua_State* L)
{
#define NC QDP_get_nc(u->links[0])
  // Check 4 or 5 arguments.
  qassert(lua_gettop(L)==4 || lua_gettop(L)==5);
  // 1. Get the gauge field.
  gauge_t *u = qopqdp_gauge_check(L, 1);
  // 2. Get time direction.
  int j_decay = luaL_checkint(L, 2);
  // 3. Get accuracy.
  double GFAccu = luaL_checknumber(L, 3);
  // 4. Get maximum number of iterations.
  int GFMax = luaL_checkint(L, 4);
	
  int OrDo = 0;
  double OrPara = 1;
	
  if (lua_gettop(L) == 5)
    {
      // 5. Overrelax?
      OrDo = 1;
      // 6. Aw man, how much do I overrelax?
      OrPara = luaL_checknumber(L, 5);
    }
	
  double tgfold;
  double tgfnew;
  double tgf_t;
  double tgf_s;
  double norm;
  int num_sdir; // Fixed because of isotropy.
  int Nd = QDP_ndim();
  int tDir = Nd - 1;
  int tdirp; // bool
  QDP_ColorMatrix *g, *new_links, *new_links2;
  g = QDP_create_M();
  new_links = QDP_create_M(); // we use this to update the link trace.
  new_links2 = QDP_create_M();
	
  // My addition:
  int mu;
  int n_gf; // Returned number of fixing terms.
	
  // Chroma has some extra logic seeing if the fixing direction
  // is the time direction because of anisotropy.
  // We don't need that. 
  if (j_decay >= 0 && j_decay < Nd) // Coulomb
    {
      if (tDir != j_decay)
	{
	  num_sdir = Nd - 2;
	  tdirp = 1;
	  norm = (double)(QDP_volume() * QLA_Nc * (num_sdir + 1));
	}
      else
	{
	  num_sdir = Nd - 1;
	  norm = (double)(QDP_volume() * QLA_Nc * num_sdir);
	  tdirp = 0;
	}
    }
  else // Landau
    {
      num_sdir = Nd - 1;
      norm = (double)(QDP_volume() * QLA_Nc * (num_sdir + 1));
      tdirp = 1;
    }
	
  // Compute initial gauge fixing term: sum(trace(U_spacelike));
  tgf_t = 0;
  tgf_s = 0;
  for (mu = 0; mu < Nd; ++mu)
    {
      if ( mu != j_decay )
	{
	  //double tgf_tmp = sum(real(trace(u[mu])));
	  QLA_ColorMatrix link_sum;
	  QDP_m_eq_sum_M(&link_sum, u->links[mu], QDP_all);
	  QLA_Real tgf_tmp;
	  QLA_R_eq_re_trace_M(&tgf_tmp, &link_sum);
	  
	  if (mu != tDir) // Keep different stock for time link trace.
	    {
	      tgf_s += tgf_tmp;
	    }
	  else
	    {
	      tgf_t += tgf_tmp;
	    }
	}
    }
	
  // Get the normalized trace sum of links.
  // We get different values depending on time dependence.
  if (tdirp) 
    {
      tgfold = (tgf_t + tgf_s)/norm;
      tgf_s = tgf_s/((double)(QDP_volume() * QLA_Nc * num_sdir));
      tgf_t = tgf_t/((double)(QDP_volume() * QLA_Nc));
    }
  else
    {
      tgf_s = tgf_s/((double)(QDP_volume() * QLA_Nc * num_sdir));
      tgfold = tgf_s;
    }
	
  // Gauge transformation matrices always start from the identity.
  QLA_Complex z;
  QLA_c_eq_r(z, 1);
  QDP_M_eq_c(g, &z, QDP_all);
	
  // Gauge fix until converged or too many iterations.
  n_gf = 0;
  int wrswitch = 1; // Switch for writing of gauge fixing term.
  double conver = 1; // Convergence criterion.
  double wreps = 0.1; // print output every reduction of conver by eps
  double cnvwr = wreps*conver; // conver for next output
	
  while ((conver > GFAccu) && n_gf < GFMax) // Loop 'til we're done!
    {
      n_gf = n_gf + 1;
		
      // Loop over checkerboards for gauge fixing.
      for (int cb = 0; cb < 2; ++cb)
	{
	  // Loop over SU(2) subgroup index.
	  for (int su2_index = 0; su2_index < QLA_Nc*(QLA_Nc-1)/2; ++su2_index)
	    {
	      // Do a gauge fixing relaxation step!
	      private_grelax(g, u, j_decay, su2_index, cb, OrDo, OrPara);
	    } 
	}
      
      // Reunitarize, borrowing from "qopqdp_gauge_makeSU("
      QDP_M_eq_funcia(g, qopqdp_makeSU, NULL, QDP_all);
      
      // Update gauge fixing term: sum(trace(U_spacelike)); 
      tgf_t = 0;
      tgf_s = 0;
      for (mu = 0; mu < Nd; ++mu)
	{
	  if ( mu != j_decay )
	    {
	      //double tgf_tmp = sum(real(trace(g(x)*u[mu](x)*g(x+mu)^dagger)));
	      // First get g(x-mu).
	      QDP_M_eq_sM(new_links, g, QDP_neighbor[mu], QDP_forward, QDP_all);
	      // Next get u[mu](x)*g(x-mu)^dagger.
	      QDP_M_eq_M_times_Ma(new_links2, u->links[mu], new_links, QDP_all);
	      // Lastly, get g(x)*u[mu](x)*g(x-mu)^dagger.
	      QDP_M_eq_M_times_M(new_links, g, new_links2, QDP_all);
	      
	      // Now get the real trace >_<
	      QLA_ColorMatrix link_sum;
	      QDP_m_eq_sum_M(&link_sum, new_links, QDP_all);
	      QLA_Real tgf_tmp;
	      QLA_R_eq_re_trace_M(&tgf_tmp, &link_sum);
	      
	      //printf0("Dir %d Trace %f\n", mu, tgf_tmp);
	      
	      if (mu != tDir) // Keep different stock for time link trace.
		{
		  tgf_s += tgf_tmp;
		}
	      else
		{
		  tgf_t += tgf_tmp;
		}
	    }
	}
      
      // Get the normalized trace sum of links.
      // We get different values depending on time dependence.
      if (tdirp) 
	{
	  tgfnew = (tgf_t + tgf_s)/norm;
	  tgf_s = tgf_s/((double)(QDP_volume() * QLA_Nc * num_sdir));
	  tgf_t = tgf_t/((double)(QDP_volume() * QLA_Nc));
	}
      else
	{
	  tgf_s = tgf_s/((double)(QDP_volume() * QLA_Nc * num_sdir));
	  tgfnew = tgf_s;
	}
      conver = fabs((tgfnew - tgfold) / tgfnew);
      
      if(wrswitch==1 && conver<cnvwr) { // Print out debugging stuff.
	printf0("COULGAUGE: iter = %d, tgfold = %f, tgfnew = %f, tgf_s = %f, tgf_t = %f\n",
		n_gf, tgfold, tgfnew, tgf_s, tgf_t);
	cnvwr *= wreps;
      }
      
      tgfold = tgfnew;
    } // End while loop
  
  // Print out final stuff.
  if ( wrswitch == 1) // Print out debugging stuff.
    {
      printf0("COULGAUGE: end: iter = %d, tgfold = %f, tgf_s = %f, tgf_t = %f\n",
	      n_gf, tgfold, tgf_s, tgf_t);
    }
  
  // Finally, gauge rotate the original matrices and overwrite.
  for (mu = 0; mu < Nd; mu++)
    {
      // First get g(x-mu).
      QDP_M_eq_sM(new_links, g, QDP_neighbor[mu], QDP_forward, QDP_all);
      // Next get u[mu](x)*g(x-mu)^dagger.
      QDP_M_eq_M_times_Ma(new_links2, u->links[mu], new_links, QDP_all);
      // Lastly, get g(x)*u[mu](x)*g(x-mu)^dagger.
      QDP_M_eq_M_times_M(u->links[mu], g, new_links2, QDP_all);
    }
  
  QDP_destroy_M(new_links);
  QDP_destroy_M(new_links2);
  QDP_destroy_M(g); // Clean it up at the end!
  
  return 1;
#undef NC
}
