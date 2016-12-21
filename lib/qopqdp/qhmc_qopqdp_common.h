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
#define qhmc_get_real_array qhmc_get_float_array
#define qhmc_push_real_array qhmc_push_float_array
#define sum_real_array QMP_sum_float_array
#define QDPOP(x) QDP_DF_ ## x
#define QDPPO(x) QDP_FD_ ## x
#else
#define qhmc_get_real_array qhmc_get_double_array
#define qhmc_push_real_array qhmc_push_double_array
#define sum_real_array QMP_sum_double_array
#define QDPOP(x) QDP_FD_ ## x
#define QDPPO(x) QDP_DF_ ## x
#endif

#if !defined(HAVE_MG) && !defined(QOP_WilsonMg)
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
  QDP_RandomState *rs;
  int nd;
  int ref;
  char *defaultPrecision;
  int defaultNc;
  int doGC;
} lattice_t;
lattice_t *qopqdp_lattice_wrap(lua_State *L, QDP_Lattice *qlat, char *defPrec, int defNc, int doGC);
lattice_t *qopqdp_create_lattice(lua_State *L, int nd, int size[], char *defPrec, int defNc);
lattice_t *qopqdp_opt_lattice(lua_State *L, int *idx, int required, lattice_t *def);
lattice_t *qopqdp_get_default_lattice(lua_State *L);
#define qopqdp_check_lattice(L,i) qopqdp_opt_lattice(L,(int[]){i},1,NULL)
#define GET_LATTICE(l) lattice_t *l = qopqdp_opt_lattice(L,&nextarg,1,NULL)
#define OPT_LATTICE(l,d) lattice_t *l = qopqdp_opt_lattice(L,&nextarg,0,d)


typedef struct {
  QDP_Subset *qgroup;
  int refcount;
} subsetGroup_t;
subsetGroup_t *qhmc_qopqdp_subsetGroup_create(lua_State *L, QDP_Subset *qgroup);
void qhmc_qopqdp_subsetGroup_free(lua_State *L, subsetGroup_t *sg);

typedef struct {
  QDP_Subset qsub;
  subsetGroup_t *group;
} subset_t;
subset_t *qopqdp_subset_create(lua_State *L, QDP_Subset qsub, subsetGroup_t *group);
subset_t *qhmc_qopqdp_opt_subset(lua_State *L, int *idx, int required, subset_t *def);
int qhmc_qopqdp_opt_as_subset_array_len(lua_State *L, int idx, int required, int def);
void qhmc_qopqdp_opt_as_subset_array(lua_State *L, int *idx, int required, int n, subset_t **t, int dn, subset_t **def);
#define qopqdp_check_lsubset(L,i) qopqdp_opt_subset(L,(int[]){i},1,NULL)
#define GET_LSUBSET(l) subset_t *l = qhmc_qopqdp_opt_subset(L,&nextarg,1,NULL)
#define OPT_LSUBSET(l,d) subset_t *l = qhmc_qopqdp_opt_subset(L,&nextarg,0,d)
QDP_Subset *qhmc_qopqdp_qsubset_from_string(lua_State *L, lattice_t *lat, const char *s, int *n);
QDP_Subset qopqdp_opt_qsubset(lua_State *L, int *idx, int reqd,
			      lattice_t *lat, QDP_Subset def);
int qopqdp_opt_as_qsubset_array_len(lua_State *L, int idx, int required, lattice_t *lat, int def);
void qopqdp_opt_as_qsubset_array(lua_State *L, int *idx, int required, lattice_t *lat, int n, QDP_Subset *t, int dn, QDP_Subset *def);
#define qopqdp_check_qsubset(L,i,l) qopqdp_opt_qsubset(L,(int[]){i},1,l,NULL)
#define GET_QSUBSET(s,l) QDP_Subset s = qopqdp_opt_qsubset(L,&nextarg,1,l,NULL)
#define OPT_QSUBSET(s,l,d) QDP_Subset s = qopqdp_opt_qsubset(L,&nextarg,0,l,d)
// n = -1 if single subset requested (not array)
#define GET_AS_QSUBSET_ARRAY(n,t,l) int n=qopqdp_opt_as_qsubset_array_len(L,nextarg,1,l,0); QDP_Subset t[abs(n)]; qopqdp_opt_as_qsubset_array(L,&nextarg,1,l,n,t,0,NULL)
#define OPT_AS_QSUBSET_ARRAY(n,t,l,dn,dt) int n=qopqdp_opt_as_qsubset_array_len(L,nextarg,0,l,dn); QDP_Subset t[abs(n)]; qopqdp_opt_as_qsubset_array(L,&nextarg,0,l,n,t,dn,dt)
#define GET_QSUBSETS(s,n,l) GET_AS_QSUBSET_ARRAY(n,s,l)
#define OPT_QSUBSETS(s,n,l,d,dn) OPT_AS_QSUBSET_ARRAY(n,s,l,dn,d)
#define GET_SUBSET GET_QSUBSET
#define OPT_SUBSET OPT_QSUBSET
#define GET_SUBSETS GET_QSUBSETS
#define OPT_SUBSETS OPT_QSUBSETS

typedef struct {
  QDP_Reader *qr;
  lattice_t *lat;
  int open;
} reader_t;
reader_t *qopqdp_reader_create(lua_State *L, const char *fn, lattice_t *lat);
reader_t *qopqdp_reader_check(lua_State *L, int idx);
void qopqdp_get_prec_type_nc(QDP_Reader *qr, int *prec, int *type, int *nc);
#define GET_READER(r) reader_t *r = qopqdp_reader_check(L,nextarg);nextarg++


typedef struct {
  QDP_Writer *qw;
  lattice_t *lat;
  int open;
  int pending;
} writer_t;
writer_t *qopqdp_writer_create(lua_State *L, const char *fn, const char *mds, lattice_t *lat);
writer_t *qopqdp_writer_check(lua_State *L, int idx);
#define GET_WRITER(r) writer_t *r = qopqdp_writer_check(L,nextarg);nextarg++


typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_RandomState *field;
  int nc;
  int doGC;
} qopqdp_rstate_t;
qopqdp_rstate_t *qopqdp_rstate_create(lua_State *L, lattice_t *lat);
qopqdp_rstate_t *qopqdp_rstate_create_unset(lua_State *L, lattice_t *lat);
qopqdp_rstate_t *qopqdp_rstate_wrap(lua_State *L, lattice_t *lat, QDP_RandomState *field, int doGC);
qopqdp_rstate_t *qopqdp_rstate_opt(lua_State *L, int *idx, int req, qopqdp_rstate_t *def);
void qhmc_qopqdp_seed_func(QDP_RandomState *r, int seed, int uniform, QDP_Subset s);
#define qopqdp_rstate_check(L,i) qopqdp_rstate_opt(L,(int[]){i},1,NULL)
#define GET_QOPQDP_RSTATE(t) qopqdp_rstate_t *t = qopqdp_rstate_opt(L,&nextarg,1,NULL)
#define OPT_QOPQDP_RSTATE(t,d) qopqdp_rstate_t *t = qopqdp_rstate_opt(L,&nextarg,0,d)
#define OPT_QOPQDP_RSTATEO(t,d) qopqdp_rstateO_t *t = qopqdp_rstateO_opt(L,&nextarg,0,d)
#define OPT_AS_QOPQDP_RSTATE_ARRAY(n,t,dn,dt) int n=qopqdp_rstate_as_array_opt_len(L,nextarg,0,dn); qopqdp_rstate_t *t[n]; qopqdp_rstate_as_array_opt(L,&nextarg,0,n,t,dn,dt)
#define OPT_QOPQDP_QRSTATE(t,d) QDP_RandomState *t=d;do{qopqdp_rstate_t *_rs=qopqdp_rstate_opt(L,&nextarg,0,NULL);if(_rs!=NULL)t=_rs->field;}while(0)
#define qopqdp_rstateO_t                 qopqdp_rstate_t
#define qopqdp_rstateO_create            qopqdp_rstate_create
#define qopqdp_rstateo_create_unset      qopqdp_rstate_create_unset
#define qopqdp_rstateO_opt               qopqdp_rstate_opt


typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_F_Real *field;
  int nc;
  int doGC;
} qopqdp_realF_t;
qopqdp_realF_t *qopqdp_realF_create(lua_State *L, lattice_t *lat);
qopqdp_realF_t *qopqdp_realF_create_unset(lua_State *L, lattice_t *lat);
qopqdp_realF_t *qopqdp_realF_opt(lua_State *L, int *idx, int req, qopqdp_realF_t *def);
typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_D_Real *field;
  int nc;
  int doGC;
} qopqdp_realD_t;
qopqdp_realD_t *qopqdp_realD_create(lua_State *L, lattice_t *lat);
qopqdp_realD_t *qopqdp_realD_create_unset(lua_State *L, lattice_t *lat);
qopqdp_realD_t *qopqdp_realD_opt(lua_State *L, int *idx, int req, qopqdp_realD_t *def);
#define qopqdp_real_check(L,i) qopqdp_real_opt(L,(int[]){i},1,NULL)
#define GET_QOPQDP_REAL(t) qopqdp_real_t *t = qopqdp_real_opt(L,&nextarg,1,NULL)
#define OPT_QOPQDP_REAL(t,d) qopqdp_real_t *t = qopqdp_real_opt(L,&nextarg,0,d)
#define OPT_QOPQDP_REALO(t,d) qopqdp_realO_t *t = qopqdp_realO_opt(L,&nextarg,0,d)
#define GET_AS_QOPQDP_REAL_ARRAY(n,t) int n=qopqdp_real_as_array_opt_len(L,nextarg,1,0); qopqdp_real_t *t[n]; qopqdp_real_as_array_opt(L,&nextarg,1,n,t,0,NULL)
#define OPT_AS_QOPQDP_REAL_ARRAY(n,t,dn,dt) int n=qopqdp_real_as_array_opt_len(L,nextarg,0,dn); qopqdp_real_t *t[n]; qopqdp_real_as_array_opt(L,&nextarg,0,n,t,dn,dt)
#if QOP_Precision == 'F'
#define qopqdp_real_t                 qopqdp_realF_t
#define qopqdp_real_create            qopqdp_realF_create
#define qopqdp_real_create_unset      qopqdp_realF_create_unset
#define qopqdp_real_wrap              qopqdp_realF_wrap
#define qopqdp_real_opt               qopqdp_realF_opt
#define qopqdp_real_as_array_opt_len  qopqdp_realF_as_array_opt_len
#define qopqdp_real_as_array_opt      qopqdp_realF_as_array_opt
#define qopqdp_realO_t                qopqdp_realD_t
#define qopqdp_realO_create           qopqdp_realD_create
#define qopqdp_realO_create_unset     qopqdp_realD_create_unset
#define qopqdp_realO_opt              qopqdp_realD_opt
#else
#define qopqdp_real_t                 qopqdp_realD_t
#define qopqdp_real_create            qopqdp_realD_create
#define qopqdp_real_create_unset      qopqdp_realD_create_unset
#define qopqdp_real_wrap              qopqdp_realD_wrap
#define qopqdp_real_opt               qopqdp_realD_opt
#define qopqdp_real_as_array_opt_len  qopqdp_realD_as_array_opt_len
#define qopqdp_real_as_array_opt      qopqdp_realD_as_array_opt
#define qopqdp_realO_t                qopqdp_realF_t
#define qopqdp_realO_create           qopqdp_realF_create
#define qopqdp_realO_create_unset     qopqdp_realF_create_unset
#define qopqdp_realO_opt              qopqdp_realF_opt
#endif


typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_F_Complex *field;
  int nc;
  int doGC;
} qopqdp_complexF_t;
qopqdp_complexF_t *qopqdp_complexF_create(lua_State *L, lattice_t *lat);
qopqdp_complexF_t *qopqdp_complexF_create_unset(lua_State *L, lattice_t *lat);
qopqdp_complexF_t *qopqdp_complexF_opt(lua_State *L, int *idx, int req, qopqdp_complexF_t *def);
typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_D_Complex *field;
  int nc;
  int doGC;
} qopqdp_complexD_t;
qopqdp_complexD_t *qopqdp_complexD_create(lua_State *L, lattice_t *lat);
qopqdp_complexD_t *qopqdp_complexD_create_unset(lua_State *L, lattice_t *lat);
qopqdp_complexD_t *qopqdp_complexD_opt(lua_State *L, int *idx, int req, qopqdp_complexD_t *def);
#define qopqdp_complex_check(L,i) qopqdp_complex_opt(L,(int[]){i},1,NULL)
#define GET_QOPQDP_COMPLEX(t) qopqdp_complex_t *t = qopqdp_complex_opt(L,&nextarg,1,NULL)
#define OPT_QOPQDP_COMPLEX(t,d) qopqdp_complex_t *t = qopqdp_complex_opt(L,&nextarg,0,d)
#define OPT_QOPQDP_COMPLEXO(t,d) qopqdp_complexO_t *t = qopqdp_complexO_opt(L,&nextarg,0,d)
#define GET_AS_QOPQDP_COMPLEX_ARRAY(n,t) int n=qopqdp_complex_as_array_opt_len(L,nextarg,1,0); qopqdp_complex_t *t[n]; qopqdp_complex_as_array_opt(L,&nextarg,1,n,t,0,NULL)
#define OPT_AS_QOPQDP_COMPLEX_ARRAY(n,t,dn,dt) int n=qopqdp_complex_as_array_opt_len(L,nextarg,0,dn); qopqdp_complex_t *t[n]; qopqdp_complex_as_array_opt(L,&nextarg,0,n,t,dn,dt)
#if QOP_Precision == 'F'
#define qopqdp_complex_t                 qopqdp_complexF_t
#define qopqdp_complex_create            qopqdp_complexF_create
#define qopqdp_complex_create_unset      qopqdp_complexF_create_unset
#define qopqdp_complex_wrap              qopqdp_complexF_wrap
#define qopqdp_complex_opt               qopqdp_complexF_opt
#define qopqdp_complex_as_array_opt_len  qopqdp_complexF_as_array_opt_len
#define qopqdp_complex_as_array_opt      qopqdp_complexF_as_array_opt
#define qopqdp_complexO_t                qopqdp_complexD_t
#define qopqdp_complexO_create           qopqdp_complexD_create
#define qopqdp_complexO_create_unset     qopqdp_complexD_create_unset
#define qopqdp_complexO_opt              qopqdp_complexD_opt
#else
#define qopqdp_complex_t                 qopqdp_complexD_t
#define qopqdp_complex_create            qopqdp_complexD_create
#define qopqdp_complex_create_unset      qopqdp_complexD_create_unset
#define qopqdp_complex_wrap              qopqdp_complexD_wrap
#define qopqdp_complex_opt               qopqdp_complexD_opt
#define qopqdp_complex_as_array_opt_len  qopqdp_complexD_as_array_opt_len
#define qopqdp_complex_as_array_opt      qopqdp_complexD_as_array_opt
#define qopqdp_complexO_t                qopqdp_complexF_t
#define qopqdp_complexO_create           qopqdp_complexF_create
#define qopqdp_complexO_create_unset     qopqdp_complexF_create_unset
#define qopqdp_complexO_opt              qopqdp_complexF_opt
#endif


typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_F_ColorVector *field;
  int nc;
  int doGC;
} qopqdp_cvectorF_t;
qopqdp_cvectorF_t *qopqdp_cvectorF_create(lua_State *L, int nc, lattice_t *lat);
qopqdp_cvectorF_t *qopqdp_cvectorF_create_unset(lua_State *L, int nc, lattice_t *lat);
qopqdp_cvectorF_t *qopqdp_cvectorF_opt(lua_State *L, int *idx, int req, qopqdp_cvectorF_t *def);
int qopqdp_cvectorF_as_array_opt_len(lua_State *L, int idx, int required, int def);
void qopqdp_cvectorF_as_array_opt(lua_State *L, int *idx, int required, int n, qopqdp_cvectorF_t **t, int dn, qopqdp_cvectorF_t **def);
typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_D_ColorVector *field;
  int nc;
  int doGC;
} qopqdp_cvectorD_t;
qopqdp_cvectorD_t *qopqdp_cvectorD_create(lua_State *L, int nc, lattice_t *lat);
qopqdp_cvectorD_t *qopqdp_cvectorD_create_unset(lua_State *L, int nc, lattice_t *lat);
qopqdp_cvectorD_t *qopqdp_cvectorD_opt(lua_State *L, int *idx, int req, qopqdp_cvectorD_t *def);
int qopqdp_cvectorD_as_array_opt_len(lua_State *L, int idx, int required, int def);
void qopqdp_cvectorD_as_array_opt(lua_State *L, int *idx, int required, int n, qopqdp_cvectorD_t **t, int dn, qopqdp_cvectorD_t **def);
#define qopqdp_cvector_check(L,i) qopqdp_cvector_opt(L,(int[]){i},1,NULL)
#define GET_QOPQDP_CVECTOR(t) qopqdp_cvector_t *t = qopqdp_cvector_opt(L,&nextarg,1,NULL)
#define OPT_QOPQDP_CVECTOR(t,d) qopqdp_cvector_t *t = qopqdp_cvector_opt(L,&nextarg,0,d)
#define OPT_QOPQDP_CVECTORO(t,d) qopqdp_cvectorO_t *t = qopqdp_cvectorO_opt(L,&nextarg,0,d)
#define GET_AS_QOPQDP_CVECTOR_ARRAY(n,t) int n=qopqdp_cvector_as_array_opt_len(L,nextarg,1,0); qopqdp_cvector_t *t[n]; qopqdp_cvector_as_array_opt(L,&nextarg,1,n,t,0,NULL)
#define OPT_AS_QOPQDP_CVECTOR_ARRAY(n,t,dn,dt) int n=qopqdp_cvector_as_array_opt_len(L,nextarg,0,dn); qopqdp_cvector_t *t[n]; qopqdp_cvector_as_array_opt(L,&nextarg,0,n,t,dn,dt)
#if QOP_Precision == 'F'
#define qopqdp_cvector_t                 qopqdp_cvectorF_t
#define qopqdp_cvector_create            qopqdp_cvectorF_create
#define qopqdp_cvector_create_unset      qopqdp_cvectorF_create_unset
#define qopqdp_cvector_wrap              qopqdp_cvectorF_wrap
#define qopqdp_cvector_opt               qopqdp_cvectorF_opt
#define qopqdp_cvector_as_array_opt_len  qopqdp_cvectorF_as_array_opt_len
#define qopqdp_cvector_as_array_opt      qopqdp_cvectorF_as_array_opt
#define qopqdp_cvectorO_t                qopqdp_cvectorD_t
#define qopqdp_cvectorO_create           qopqdp_cvectorD_create
#define qopqdp_cvectorO_create_unset     qopqdp_cvectorD_create_unset
#define qopqdp_cvectorO_opt              qopqdp_cvectorD_opt
#else
#define qopqdp_cvector_t                 qopqdp_cvectorD_t
#define qopqdp_cvector_create            qopqdp_cvectorD_create
#define qopqdp_cvector_create_unset      qopqdp_cvectorD_create_unset
#define qopqdp_cvector_wrap              qopqdp_cvectorD_wrap
#define qopqdp_cvector_opt               qopqdp_cvectorD_opt
#define qopqdp_cvector_as_array_opt_len  qopqdp_cvectorD_as_array_opt_len
#define qopqdp_cvector_as_array_opt      qopqdp_cvectorD_as_array_opt
#define qopqdp_cvectorO_t                qopqdp_cvectorF_t
#define qopqdp_cvectorO_create           qopqdp_cvectorF_create
#define qopqdp_cvectorO_create_unset     qopqdp_cvectorF_create_unset
#define qopqdp_cvectorO_opt              qopqdp_cvectorF_opt
#endif


typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_F_DiracFermion *field;
  int nc;
  int doGC;
} qopqdp_dfermionF_t;
qopqdp_dfermionF_t *qopqdp_dfermionF_create(lua_State *L, int nc, lattice_t *lat);
qopqdp_dfermionF_t *qopqdp_dfermionF_create_unset(lua_State *L, int nc, lattice_t *lat);
qopqdp_dfermionF_t *qopqdp_dfermionF_opt(lua_State *L, int *idx, int req, qopqdp_dfermionF_t *def);
int qopqdp_dfermionF_as_array_opt_len(lua_State *L, int idx, int required, int def);
typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_D_DiracFermion *field;
  int nc;
  int doGC;
} qopqdp_dfermionD_t;
qopqdp_dfermionD_t *qopqdp_dfermionD_create(lua_State *L, int nc, lattice_t *lat);
qopqdp_dfermionD_t *qopqdp_dfermionD_create_unset(lua_State *L, int nc, lattice_t *lat);
qopqdp_dfermionD_t *qopqdp_dfermionD_opt(lua_State *L, int *idx, int req, qopqdp_dfermionD_t *def);
int qopqdp_dfermionD_as_array_opt_len(lua_State *L, int idx, int required, int def);
void qopqdp_dfermionD_as_array_opt(lua_State *L, int *idx, int required, int n, qopqdp_dfermionD_t **t, int dn, qopqdp_dfermionD_t **def);
#define qopqdp_dfermion_check(L,i) qopqdp_dfermion_opt(L,(int[]){i},1,NULL)
#define GET_QOPQDP_DFERMION(t) qopqdp_dfermion_t *t = qopqdp_dfermion_opt(L,&nextarg,1,NULL)
#define OPT_QOPQDP_DFERMION(t,d) qopqdp_dfermion_t *t = qopqdp_dfermion_opt(L,&nextarg,0,d)
#define OPT_QOPQDP_DFERMIONO(t,d) qopqdp_dfermionO_t *t = qopqdp_dfermionO_opt(L,&nextarg,0,d)
#define GET_AS_QOPQDP_DFERMION_ARRAY(n,t) int n=qopqdp_dfermion_as_array_opt_len(L,nextarg,1,0); qopqdp_dfermion_t *t[n]; qopqdp_dfermion_as_array_opt(L,&nextarg,1,n,t,0,NULL)
#define OPT_AS_QOPQDP_DFERMION_ARRAY(n,t,dn,dt) int n=qopqdp_dfermion_as_array_opt_len(L,nextarg,0,dn); qopqdp_dfermion_t *t[n]; qopqdp_dfermion_as_array_opt(L,&nextarg,0,n,t,dn,dt)
#if QOP_Precision == 'F'
#define qopqdp_dfermion_t                 qopqdp_dfermionF_t
#define qopqdp_dfermion_create            qopqdp_dfermionF_create
#define qopqdp_dfermion_create_unset      qopqdp_dfermionF_create_unset
#define qopqdp_dfermion_wrap              qopqdp_dfermionF_wrap
#define qopqdp_dfermion_opt               qopqdp_dfermionF_opt
#define qopqdp_dfermion_as_array_opt_len  qopqdp_dfermionF_as_array_opt_len
#define qopqdp_dfermion_as_array_opt      qopqdp_dfermionF_as_array_opt
#define qopqdp_dfermionO_t                qopqdp_dfermionD_t
#define qopqdp_dfermionO_create           qopqdp_dfermionD_create
#define qopqdp_dfermionO_create_unset     qopqdp_dfermionD_create_unset
#define qopqdp_dfermionO_opt              qopqdp_dfermionD_opt
#else
#define qopqdp_dfermion_t                 qopqdp_dfermionD_t
#define qopqdp_dfermion_create            qopqdp_dfermionD_create
#define qopqdp_dfermion_create_unset      qopqdp_dfermionD_create_unset
#define qopqdp_dfermion_wrap              qopqdp_dfermionD_wrap
#define qopqdp_dfermion_opt               qopqdp_dfermionD_opt
#define qopqdp_dfermion_as_array_opt_len  qopqdp_dfermionD_as_array_opt_len
#define qopqdp_dfermion_as_array_opt      qopqdp_dfermionD_as_array_opt
#define qopqdp_dfermionO_t                qopqdp_dfermionF_t
#define qopqdp_dfermionO_create           qopqdp_dfermionF_create
#define qopqdp_dfermionO_create_unset     qopqdp_dfermionF_create_unset
#define qopqdp_dfermionO_opt              qopqdp_dfermionF_opt
#endif


typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_F_ColorMatrix *field;
  int nc;
  int doGC;
} qopqdp_cmatrixF_t;
qopqdp_cmatrixF_t *qopqdp_cmatrixF_wrap(lua_State *L, lattice_t *lat, QDP_F_ColorMatrix *field, int doGC);
qopqdp_cmatrixF_t *qopqdp_cmatrixF_create(lua_State *L, int nc, lattice_t *lat);
qopqdp_cmatrixF_t *qopqdp_cmatrixF_create_unset(lua_State *L, int nc, lattice_t *lat);
qopqdp_cmatrixF_t *qopqdp_cmatrixF_opt(lua_State *L, int *idx, int req, qopqdp_cmatrixF_t *def);
typedef struct {
  lattice_t *lat;
  QDP_Lattice *qlat;
  QDP_D_ColorMatrix *field;
  int nc;
  int doGC;
} qopqdp_cmatrixD_t;
qopqdp_cmatrixD_t *qopqdp_cmatrixD_wrap(lua_State *L, lattice_t *lat, QDP_D_ColorMatrix *field, int doGC);
qopqdp_cmatrixD_t *qopqdp_cmatrixD_create(lua_State *L, int nc, lattice_t *lat);
qopqdp_cmatrixD_t *qopqdp_cmatrixD_create_unset(lua_State *L, int nc, lattice_t *lat);
qopqdp_cmatrixD_t *qopqdp_cmatrixD_opt(lua_State *L, int *idx, int req, qopqdp_cmatrixD_t *def);
#define qopqdp_cmatrix_check(L,i) qopqdp_cmatrix_opt(L,(int[]){i},1,NULL)
#define GET_QOPQDP_CMATRIX(t) qopqdp_cmatrix_t *t = qopqdp_cmatrix_opt(L,&nextarg,1,NULL)
#define OPT_QOPQDP_CMATRIX(t,d) qopqdp_cmatrix_t *t = qopqdp_cmatrix_opt(L,&nextarg,0,d)
#define OPT_QOPQDP_CMATRIXO(t,d) qopqdp_cmatrixO_t *t = qopqdp_cmatrixO_opt(L,&nextarg,0,d)
#define GET_AS_QOPQDP_CMATRIX_ARRAY(n,t) int n=qopqdp_cmatrix_as_array_opt_len(L,nextarg,1,0); qopqdp_cmatrix_t *t[n]; qopqdp_cmatrix_as_array_opt(L,&nextarg,1,n,t,0,NULL)
#define OPT_AS_QOPQDP_CMATRIX_ARRAY(n,t,dn,dt) int n=qopqdp_cmatrix_as_array_opt_len(L,nextarg,0,dn); qopqdp_cmatrix_t *t[n]; qopqdp_cmatrix_as_array_opt(L,&nextarg,0,n,t,dn,dt)
#if QOP_Precision == 'F'
#define qopqdp_cmatrix_t                 qopqdp_cmatrixF_t
#define qopqdp_cmatrix_create            qopqdp_cmatrixF_create
#define qopqdp_cmatrix_create_unset      qopqdp_cmatrixF_create_unset
#define qopqdp_cmatrix_wrap              qopqdp_cmatrixF_wrap
#define qopqdp_cmatrix_opt               qopqdp_cmatrixF_opt
#define qopqdp_cmatrix_as_array_opt_len  qopqdp_cmatrixF_as_array_opt_len
#define qopqdp_cmatrix_as_array_opt      qopqdp_cmatrixF_as_array_opt
#define qopqdp_cmatrixO_t                qopqdp_cmatrixD_t
#define qopqdp_cmatrixO_create           qopqdp_cmatrixD_create
#define qopqdp_cmatrixO_create_unset     qopqdp_cmatrixD_create_unset
#define qopqdp_cmatrixO_opt              qopqdp_cmatrixD_opt
#else
#define qopqdp_cmatrix_t                 qopqdp_cmatrixD_t
#define qopqdp_cmatrix_create            qopqdp_cmatrixD_create
#define qopqdp_cmatrix_create_unset      qopqdp_cmatrixD_create_unset
#define qopqdp_cmatrix_wrap              qopqdp_cmatrixD_wrap
#define qopqdp_cmatrix_opt               qopqdp_cmatrixD_opt
#define qopqdp_cmatrix_as_array_opt_len  qopqdp_cmatrixD_as_array_opt_len
#define qopqdp_cmatrix_as_array_opt      qopqdp_cmatrixD_as_array_opt
#define qopqdp_cmatrixO_t                qopqdp_cmatrixF_t
#define qopqdp_cmatrixO_create           qopqdp_cmatrixF_create
#define qopqdp_cmatrixO_create_unset     qopqdp_cmatrixF_create_unset
#define qopqdp_cmatrixO_opt              qopqdp_cmatrixF_opt
#endif


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
#define GROUP_TAH (GROUP_T+GROUP_AH)
#define GROUP_TYPE 3 // (g&GROUP_TYPE)->{GL,U,H,AH}
typedef struct {
  double time;
  double flops;
  lattice_t *lat;
  QDP_Lattice *qlat;
  int nd, nc;
  QDP_F_ColorMatrix *links[];
} gaugeF_t;
gaugeF_t *qopqdp_gaugeF_create(lua_State *L, int nc, lattice_t *lat);
gaugeF_t *qopqdp_gaugeF_opt(lua_State *L, int *idx,int required,gaugeF_t *def);
gaugeF_t *qopqdp_gaugeF_check(lua_State *L, int idx);
void qopqdp_gaugeF_array_check(lua_State *L, int idx, int n, gaugeF_t *g[n]);
int qopqdp_gaugeF_coulomb(lua_State *L);
void qopqdp_F_makeSU(NCPROT QLA_F_ColorMatrix(*m), int idx, void *args);
typedef struct {
  double time;
  double flops;
  lattice_t *lat;
  QDP_Lattice *qlat;
  int nd, nc, group;
  QDP_D_ColorMatrix *links[];
} gaugeD_t;
gaugeD_t *qopqdp_gaugeD_create(lua_State *L, int nc, lattice_t *lat);
gaugeD_t *qopqdp_gaugeD_opt(lua_State *L, int *idx,int required,gaugeD_t *def);
gaugeD_t *qopqdp_gaugeD_check(lua_State *L, int idx);
void qopqdp_gaugeD_array_check(lua_State *L, int idx, int n, gaugeD_t *g[n]);
int qopqdp_gaugeD_coulomb(lua_State *L);
void qopqdp_D_makeSU(NCPROT QLA_D_ColorMatrix(*m), int idx, void *args);
#if QOP_Precision == 'F'
#define gauge_t gaugeF_t
#define gaugeO_t gaugeD_t
#define qopqdp_gauge_create qopqdp_gaugeF_create
#define qopqdp_gaugeO_create qopqdp_gaugeD_create
#define qopqdp_gauge_opt qopqdp_gaugeF_opt
#define qopqdp_gauge_check qopqdp_gaugeF_check
#define qopqdp_gauge_array_check qopqdp_gaugeF_array_check
#define qopqdp_gaugeO_array_check qopqdp_gaugeD_array_check
#define qopqdp_gauge_coulomb qopqdp_gaugeF_coulomb
#define qopqdp_makeSU qopqdp_F_makeSU
#else
#define gauge_t gaugeD_t
#define gaugeO_t gaugeF_t
#define qopqdp_gauge_create qopqdp_gaugeD_create
#define qopqdp_gaugeO_create qopqdp_gaugeF_create
#define qopqdp_gauge_opt qopqdp_gaugeD_opt
#define qopqdp_gauge_check qopqdp_gaugeD_check
#define qopqdp_gauge_array_check qopqdp_gaugeD_array_check
#define qopqdp_gaugeO_array_check qopqdp_gaugeF_array_check
#define qopqdp_gauge_coulomb qopqdp_gaugeD_coulomb
#define qopqdp_makeSU qopqdp_D_makeSU
#endif
#define GET_GAUGE(g) gauge_t *g = qopqdp_gauge_check(L,nextarg++)
#define OPT_GAUGE(g,d) gauge_t *g = qopqdp_gauge_opt(L,&nextarg,0,d)
#define GET_GAUGE_COEFFS(c) QOP_gauge_coeffs_t c; get_gauge_coeffs(L,&c,nextarg++)


typedef struct {
  double time;
  double flops;
  double rsq;
  int its;
  int *r0;
  int nc;
  lattice_t *lat;
  QOP_bc_t bc;
  QOP_staggered_sign_t ssign;
  QOP_asqtad_coeffs_t coeffs;
  QOP_FermionLinksAsqtad *fl;
  QOP_F_FermionLinksAsqtad *ffl;
  gauge_t *g;
  gauge_t *gLong;
  QDP_Complex **fatphase, **longphase;
} asqtad_t;
asqtad_t *qopqdp_asqtad_create(lua_State *L, int nc, lattice_t *lat);
asqtad_t *qopqdp_asqtad_check(lua_State *L, int idx);
#define GET_ASQTAD(a) asqtad_t *a = qopqdp_asqtad_check(L,nextarg++)


typedef struct {
  double time;
  double flops;
  double rsq;
  int its;
  int nc;
  lattice_t *lat;
  QOP_hisq_coeffs_t coeffs;
  QOP_FermionLinksHisq *fl;
  QOP_F_FermionLinksHisq *ffl;
  gauge_t *g;
  double u0;
  double f7lf;
  QDP_Complex **fatphase, **longphase;
} hisq_t;
hisq_t *qopqdp_hisq_create(lua_State *L, int nc, lattice_t *lat);
hisq_t *qopqdp_hisq_check(lua_State *L, int idx);
#define GET_HISQ(a) hisq_t *a = qopqdp_hisq_check(L,nextarg++)


void asqtadInvert(QOP_info_t *info, QOP_FermionLinksAsqtad *fla,
		  QOP_invert_arg_t *invarg, QOP_resid_arg_t *residarg[],
		  QLA_Real *masses, int nm, QDP_ColorVector *prop[],
		  QDP_ColorVector *source);


typedef struct {
  double time;
  double flops;
  double rsq;
  int its;
  int nc;
  lattice_t *lat;
  QOP_wilson_coeffs_t coeffs;
  QOP_FermionLinksWilson *fl;
  QOP_F_FermionLinksWilson *ffl;
  QOP_WilsonMg *mg;
  gauge_t *g;
} wilson_t;
wilson_t *qopqdp_wilson_create(lua_State *L, int nc, lattice_t *lat);
wilson_t *qopqdp_wilson_check(lua_State *L, int idx);
#define GET_WILSON(w) wilson_t *w = qopqdp_wilson_check(L,nextarg);nextarg++


typedef struct {
  double time;
  double flops;
  double rsq;
  int its;
  int ls;
  QOP_dw_coeffs_t coeffs;
  QOP_FermionLinksDW *fl;
  QOP_F_FermionLinksDW *ffl;
  gauge_t *g;
} dw_t;
dw_t *qopqdp_dw_create(lua_State *L, int nc, lattice_t *lat);
dw_t *qopqdp_dw_check(lua_State *L, int idx);
#define GET_DW(d) dw_t *d = qopqdp_dw_check(L,nextarg++)


QDP_Subset *qhmcqdp_get_timeslices(lattice_t *lat);
QDP_Subset *qhmcqdp_get_eodir(lattice_t *lat, int dir);
QOP_evenodd_t qopqdp_opt_evenodd(lua_State *L, int *idx, int required, QOP_evenodd_t def);
#define qopqdp_check_evenodd(L,i) qopqdp_opt_evenodd(L,(int[]){i},1,QOP_EVENODD)
#define OPT_EVENODD(t,d) QOP_evenodd_t t = qopqdp_opt_evenodd(L,&nextarg,0,d)

void qhmc_qopqdp_getCopyHyper(QDP_Shift *map, QDP_Subset **subset,
			      QDP_Lattice *rlat, int roff[], int rlen[], int sdir[],
			      QDP_Lattice *slat, int soff[], int num);

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

void open_qopqdp_smear(lua_State *L);
