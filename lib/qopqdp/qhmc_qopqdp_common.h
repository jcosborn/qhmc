#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifndef QOP_Precision
#define QOP_Precision 'D'
#define QDP_Precision 'D'
#define QLA_Precision 'D'
#endif
#include <qop.h>
#include "qhmc_config_internal.h"
#ifdef HAVE_NC3
#include <qop_f3.h>
#include <qop_d3.h>
#include <qop_df3.h>
#endif
#include "qhmc_internal.h"

#if QOP_Colors == 'N'

//#define NC nc
#undef QOP_Nc
#define QOP_Nc NC
#undef QDP_Nc
#define QDP_Nc QOP_Nc
#undef QLA_Nc
#define QLA_Nc QOP_Nc
#define NCPROT int NC,
#define NCPROTVOID int NC
#define NCARG NC,
#define NCARGVOID NC

#define PARENIFNOTEMPTY(x) IFEMPTY(x,(x),x)
#define IFEMPTY(t,f,x) IFEMPTY_(COMMAFUNC x (), t, f, dum)
#define IFEMPTY_(...) THIRDARG(__VA_ARGS__)
#define THIRDARG(a,b,c,...) c
#define COMMAFUNC(...) ,

#ifndef QLA_ColorMatrix
#if QOP_Precision == 'F'
#define QLA_ColorVector(x)  QLA_FN_ColorVector (QLA_Nc,PARENIFNOTEMPTY(x))
#define QLA_ColorMatrix(x)  QLA_FN_ColorMatrix (QLA_Nc,PARENIFNOTEMPTY(x))
#define QLA_DiracFermion(x) QLA_FN_DiracFermion(QLA_Nc,PARENIFNOTEMPTY(x))
#define QLA_N_ColorVector(n,x)  QLA_FN_ColorVector (n,x)
#define QLA_N_ColorMatrix(n,x)  QLA_FN_ColorMatrix (n,x)
#define QLA_N_DiracFermion(n,x) QLA_FN_DiracFermion(n,x)
#else
#define QLA_ColorVector(x)  QLA_DN_ColorVector (QLA_Nc,PARENIFNOTEMPTY(x))
#define QLA_ColorMatrix(x)  QLA_DN_ColorMatrix (QLA_Nc,PARENIFNOTEMPTY(x))
#define QLA_DiracFermion(x) QLA_DN_DiracFermion(QLA_Nc,PARENIFNOTEMPTY(x))
#define QLA_N_ColorVector(n,x)  QLA_DN_ColorVector (n,x)
#define QLA_N_ColorMatrix(n,x)  QLA_DN_ColorMatrix (n,x)
#define QLA_N_DiracFermion(n,x) QLA_DN_DiracFermion(n,x)
#endif
#define QLA_F_ColorVector(x)  QLA_FN_ColorVector (QLA_Nc,PARENIFNOTEMPTY(x))
#define QLA_F_ColorMatrix(x)  QLA_FN_ColorMatrix (QLA_Nc,PARENIFNOTEMPTY(x))
#define QLA_F_DiracFermion(x) QLA_FN_DiracFermion(QLA_Nc,PARENIFNOTEMPTY(x))
#define QLA_D_ColorVector(x)  QLA_DN_ColorVector (QLA_Nc,PARENIFNOTEMPTY(x))
#define QLA_D_ColorMatrix(x)  QLA_DN_ColorMatrix (QLA_Nc,PARENIFNOTEMPTY(x))
#define QLA_D_DiracFermion(x) QLA_DN_DiracFermion(QLA_Nc,PARENIFNOTEMPTY(x))
#endif

#else //#if QOP_Colors == 'N'

#define NCPROT
#define NCPROTVOID void
#define NCARG
#define NCARGVOID
#ifndef QLA_ColorMatrix
#define QLA_ColorMatrix(x) QLA_ColorMatrix x
#define QLA_F_ColorMatrix(x) QLA_F_ColorMatrix x
#define QLA_D_ColorMatrix(x) QLA_D_ColorMatrix x
#endif

#endif //#if QOP_Colors == 'N'


#ifndef QLA_N_ColorMatrix
#if QOP_Precision == 'F'
#define QLA_N_ColorVector(n,x)  QLA_FN_ColorVector (n,x)
#define QLA_N_ColorMatrix(n,x)  QLA_FN_ColorMatrix (n,x)
#define QLA_N_DiracFermion(n,x) QLA_FN_DiracFermion(n,x)
#else
#define QLA_N_ColorVector(n,x)  QLA_DN_ColorVector (n,x)
#define QLA_N_ColorMatrix(n,x)  QLA_DN_ColorMatrix (n,x)
#define QLA_N_DiracFermion(n,x) QLA_DN_DiracFermion(n,x)
#endif
#endif

#if QOP_Precision == 'F'
#define get_real_array get_float_array
#define push_real_array push_float_array
#define QDPOP(x) QDP_DF_ ## x
#define QDPPO(x) QDP_FD_ ## x
#else
#define get_real_array get_double_array
#define push_real_array push_double_array
#define QDPOP(x) QDP_FD_ ## x
#define QDPPO(x) QDP_DF_ ## x
#endif

#ifndef HAVE_MG
#define QOP_WilsonMg void
#define QOP_wilsonMgNew(...) NULL
#define QOP_wilsonMgSet(...) ((void)0)
#define QOP_wilsonMgSetArray(...) ((void)0)
#define QOP_wilsonMgSetLinks(...) ((void)0)
#define QOP_wilsonMgSetup(...) ((void)0)
#define QOP_D3_wilsonMgSolve(...) ((void)0)
#endif


#undef if0
#define if0 if(QDP_this_node==0)

extern QLA_RandomState qopqdp_nrs;
extern QDP_RandomState *qopqdp_srs;


typedef struct {
  QDP_Lattice *qlat;
  QDP_Subset *timeslices;
  QDP_Subset *staggered;
  QDP_Subset **eodir;
  int nd;
  int ref;
} lattice_t;
lattice_t *qopqdp_lattice_create(lua_State *L, int nd, int size[]);
lattice_t *qopqdp_opt_lattice(lua_State *L, int *idx, int required, lattice_t *def);
lattice_t *qopqdp_get_default_lattice(lua_State *L);
#define qopqdp_check_lattice(L,i) qopqdp_opt_lattice(L,(int[]){i},1,NULL)
#define GET_LATTICE(l) lattice_t *l = qopqdp_opt_lattice(L,&nextarg,1,NULL)
#define OPT_LATTICE(l,d) lattice_t *l = qopqdp_opt_lattice(L,&nextarg,0,d)


typedef struct {
  QDP_Reader *qr;
} reader_t;
reader_t *qopqdp_reader_create(lua_State* L, const char *fn, lattice_t *lat);
reader_t *qopqdp_reader_check(lua_State *L, int idx);
void qopqdp_get_prec_type_nc(QDP_Reader *qr, int *prec, int *type, int *nc);


typedef struct {
  QDP_Writer *qw;
} writer_t;
writer_t *qopqdp_writer_create(lua_State* L, const char *fn, const char *mds);
writer_t *qopqdp_writer_check(lua_State *L, int idx);


typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_F_Complex *c;
} cscalarF_t;
cscalarF_t *qopqdp_cscalarF_create(lua_State* L, lattice_t *lat);
cscalarF_t *qopqdp_cscalarF_check(lua_State *L, int idx);
void qopqdp_cscalarF_array_check(lua_State *L, int idx, int n, cscalarF_t *s[n]);
typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_D_Complex *c;
} cscalarD_t;
cscalarD_t *qopqdp_cscalarD_create(lua_State* L, lattice_t *lat);
cscalarD_t *qopqdp_cscalarD_check(lua_State *L, int idx);
void qopqdp_cscalarD_array_check(lua_State *L, int idx, int n, cscalarD_t *s[n]);
#define GET_CSCALAR(s) cscalar_t *s = qopqdp_cscalar_check(L,nextarg); nextarg++
#if QOP_Precision == 'F'
#define cscalar_t cscalarF_t
#define qopqdp_cscalar_create qopqdp_cscalarF_create
#define qopqdp_cscalar_check qopqdp_cscalarF_check
#define qopqdp_cscalar_array_check qopqdp_cscalarF_array_check
#else
#define cscalar_t cscalarD_t
#define qopqdp_cscalar_create qopqdp_cscalarD_create
#define qopqdp_cscalar_check qopqdp_cscalarD_check
#define qopqdp_cscalar_array_check qopqdp_cscalarD_array_check
#endif


#if 0 // not currently used
#define GROUP_GL   0 // general
#define GROUP_U    1 // unitary
#define GROUP_H    2 // Hermitian
#define GROUP_AH   3 // anti-Hermitian
#define GROUP_S    4 // special
#define GROUP_T    8 // traceless
#define GROUP_SL  (GROUP_S+GROUP_GL)
#define GROUP_SU  (GROUP_S+GROUP_U)
#define GROUP_TGL (GROUP_T+GROUP_GL)
#define GROUP_TH  (GROUP_T+GROUP_H)
#define GROUP_TAH (GROUP_T+GROUP_H)
#endif
typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  int nd, nc;
  QDP_F_ColorMatrix *links[];
} gaugeF_t;
gaugeF_t *qopqdp_gaugeF_create(lua_State *L, int nc, lattice_t *lat);
gaugeF_t *qopqdp_gaugeF_check(lua_State *L, int idx);
void qopqdp_gaugeF_array_check(lua_State *L, int idx, int n, gaugeF_t *g[n]);
int qopqdp_gaugeF_coulomb(lua_State* L);
void qopqdp_F_makeSU(NCPROT QLA_F_ColorMatrix(*m), int idx, void *args);
typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  int nd, nc, group;
  QDP_D_ColorMatrix *links[];
} gaugeD_t;
gaugeD_t *qopqdp_gaugeD_create(lua_State *L, int nc, lattice_t *lat);
gaugeD_t *qopqdp_gaugeD_check(lua_State *L, int idx);
void qopqdp_gaugeD_array_check(lua_State *L, int idx, int n, gaugeD_t *g[n]);
int qopqdp_gaugeD_coulomb(lua_State* L);
void qopqdp_D_makeSU(NCPROT QLA_D_ColorMatrix(*m), int idx, void *args);
#if QOP_Precision == 'F'
#define gauge_t gaugeF_t
#define gaugeO_t gaugeD_t
#define qopqdp_gauge_create qopqdp_gaugeF_create
#define qopqdp_gaugeO_create qopqdp_gaugeD_create
#define qopqdp_gauge_check qopqdp_gaugeF_check
#define qopqdp_gauge_array_check qopqdp_gaugeF_array_check
#define qopqdp_gauge_coulomb qopqdp_gaugeF_coulomb
#define qopqdp_makeSU qopqdp_F_makeSU
#else
#define gauge_t gaugeD_t
#define gaugeO_t gaugeF_t
#define qopqdp_gauge_create qopqdp_gaugeD_create
#define qopqdp_gaugeO_create qopqdp_gaugeF_create
#define qopqdp_gauge_check qopqdp_gaugeD_check
#define qopqdp_gauge_array_check qopqdp_gaugeD_array_check
#define qopqdp_gauge_coulomb qopqdp_gaugeD_coulomb
#define qopqdp_makeSU qopqdp_D_makeSU
#endif
#define GET_GAUGE(g) gauge_t *g = qopqdp_gauge_check(L,nextarg++)
#define GET_GAUGE_COEFFS(c) QOP_gauge_coeffs_t c; get_gauge_coeffs(L,&c,nextarg++)


typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  double time;
  double flops;
  int nd, nc;
  QDP_F_ColorMatrix *force[];
} forceF_t;
forceF_t *qopqdp_forceF_create(lua_State *L, int nc, lattice_t *lat);
forceF_t *qopqdp_forceF_check(lua_State *L, int idx);
void qopqdp_forceF_array_check(lua_State *L, int idx, int n, forceF_t *g[n]);
typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  double time;
  double flops;
  int nd, nc;
  QDP_D_ColorMatrix *force[];
} forceD_t;
forceD_t *qopqdp_forceD_create(lua_State *L, int nc, lattice_t *lat);
forceD_t *qopqdp_forceD_check(lua_State *L, int idx);
void qopqdp_forceD_array_check(lua_State *L, int idx, int n, forceD_t *g[n]);
#if QOP_Precision == 'F'
#define force_t forceF_t
#define qopqdp_force_create qopqdp_forceF_create
#define qopqdp_force_check qopqdp_forceF_check
#define qopqdp_force_array_check qopqdp_forceF_array_check
#else
#define force_t forceD_t
#define qopqdp_force_create qopqdp_forceD_create
#define qopqdp_force_check qopqdp_forceD_check
#define qopqdp_force_array_check qopqdp_forceD_array_check
#endif


typedef struct {
  double time;
  double flops;
  int its;
  int nd;
  int *r0;
  QOP_bc_t bc;
  QOP_staggered_sign_t ssign;
  QOP_asqtad_coeffs_t coeffs;
  QOP_FermionLinksAsqtad *fl;
  QOP_F_FermionLinksAsqtad *ffl;
  gauge_t *g;
  QDP_Complex **fatphase, **longphase;
} asqtad_t;
asqtad_t *qopqdp_asqtad_create(lua_State *L);
asqtad_t *qopqdp_asqtad_check(lua_State *L, int idx);


typedef struct {
  double time;
  double flops;
  int its;
  QOP_hisq_coeffs_t coeffs;
  QOP_FermionLinksHisq *fl;
  QOP_F_FermionLinksHisq *ffl;
  gauge_t *g;
  double u0;
  double f7lf;
  QDP_Complex **fatphase, **longphase;
} hisq_t;
hisq_t *qopqdp_hisq_create(lua_State *L);
hisq_t *qopqdp_hisq_check(lua_State *L, int idx);


typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_ColorVector *cv;
  int nc;
} squark_t;
squark_t *qopqdp_squark_create(lua_State *L, int nc, lattice_t *lat);
squark_t *qopqdp_squark_create_unset(lua_State *L, int nc, lattice_t *lat);
squark_t *qopqdp_squark_check(lua_State *L, int idx);
void qopqdp_squark_array_check(lua_State *L, int idx, int n, squark_t *q[n]);
#define GET_SQUARK(s) squark_t *s = qopqdp_squark_check(L,nextarg++)

void asqtadInvert(QOP_info_t *info, QOP_FermionLinksAsqtad *fla,
		  QOP_invert_arg_t *invarg, QOP_resid_arg_t *residarg[],
		  QLA_Real *masses, int nm, QDP_ColorVector *prop[],
		  QDP_ColorVector *source);


typedef struct {
  double time;
  double flops;
  int its;
  QOP_wilson_coeffs_t coeffs;
  QOP_FermionLinksWilson *fl;
  QOP_F_FermionLinksWilson *ffl;
  QOP_WilsonMg *mg;
  gauge_t *g;
} wilson_t;
wilson_t *qopqdp_wilson_create(lua_State *L);
wilson_t *qopqdp_wilson_check(lua_State *L, int idx);


typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_DiracFermion *df;
  int nc;
} wquark_t;
wquark_t *qopqdp_wquark_create(lua_State *L, int nc, lattice_t *lat);
wquark_t *qopqdp_wquark_create_unset(lua_State *L, int nc, lattice_t *lat);
wquark_t *qopqdp_wquark_check(lua_State *L, int idx);
void qopqdp_wquark_array_check(lua_State *L, int idx, int n, wquark_t *q[n]);
#define GET_WQUARK(w) wquark_t *w = qopqdp_wquark_check(L,nextarg++)


typedef struct {
  double time;
  double flops;
  int its;
  int ls;
  QOP_dw_coeffs_t coeffs;
  QOP_FermionLinksDW *fl;
  QOP_F_FermionLinksDW *ffl;
  gauge_t *g;
} dw_t;
dw_t *qopqdp_dw_create(lua_State *L);
dw_t *qopqdp_dw_check(lua_State *L, int idx);


typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_DiracFermion **df;
  int ls, nc;
} dwquark_t;
dwquark_t *qopqdp_dwquark_create(lua_State *L, int ls, int nc, lattice_t *lat);
dwquark_t *qopqdp_dwquark_check(lua_State *L, int idx);
void qopqdp_dwquark_array_check(lua_State *L, int idx, int n, dwquark_t *q[n]);
#define GET_DWQUARK(d) dwquark_t *d = qopqdp_dwquark_check(L,nextarg++)


QDP_Subset *qhmcqdp_get_timeslices(lattice_t *lat);
QDP_Subset *qhmcqdp_get_eodir(lattice_t *lat, int dir);
QOP_evenodd_t qopqdp_check_evenodd(lua_State *L, int idx);
QDP_Subset qopqdp_opt_subset(lua_State *L, int *idx, int reqd,
			     lattice_t *lat, QDP_Subset def);
QDP_Subset *qopqdp_opt_subsets(lua_State *L, int *idx, int reqd,
			       lattice_t *lat, QDP_Subset def[], int *nsub);
#define qopqdp_check_subset(L,i,l) qopqdp_opt_subset(L,(int[]){i},1,l,NULL)
#define GET_SUBSET(s,l) QDP_Subset s = qopqdp_opt_subset(L,&nextarg,1,l,NULL)
#define OPT_SUBSET(s,l,d) QDP_Subset s = qopqdp_opt_subset(L,&nextarg,0,l,d)
#define GET_SUBSETS(s,n,l) int n=1; QDP_Subset *s = qopqdp_opt_subsets(L,&nextarg,1,l,NULL,&n)
#define OPT_SUBSETS(s,n,l,d,dn) int n=dn; QDP_Subset *s = qopqdp_opt_subsets(L,&nextarg,0,l,d,&n)

QLA_Real infnorm_M(QDP_ColorMatrix *m, QDP_Subset s);
int check_uniform_M(QDP_ColorMatrix *m, QDP_Subset s);
void sqrt_deriv(QDP_ColorMatrix *deriv, QDP_ColorMatrix *sqrtM,
		QDP_ColorMatrix *M, QDP_ColorMatrix *chain, QDP_Subset sub);
void projectU_deriv(QDP_ColorMatrix *deriv, QDP_ColorMatrix *proj,
		    QDP_ColorMatrix *mat, QDP_ColorMatrix *chain, QDP_Subset sub);
void sylsolve_site(NCPROT QLA_ColorMatrix(*x), QLA_ColorMatrix(*a),
		   QLA_ColorMatrix(*b), QLA_ColorMatrix(*c));
void exp_deriv_site(NCPROT QLA_ColorMatrix *deriv, QLA_Real *r,
		    QLA_ColorMatrix *M, QLA_ColorMatrix *chain);
void exp_deriv(QDP_ColorMatrix *deriv, QLA_Real *r,
	       QDP_ColorMatrix *M, QDP_ColorMatrix *chain, QDP_Subset sub);

void open_qopqdp_smear(lua_State* L);
