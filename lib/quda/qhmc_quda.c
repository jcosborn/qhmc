#include "qhmc_internal.h"
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
//#include <qmp.h>
#include <quda.h>

/*
typedef int (*QudaCommsMap)(const int *coords, void *fdata);
void initCommsGridQuda(int nDim, const int *dims, QudaCommsMap func, void *fdata);
QudaInvertParam newQudaInvertParam(void);
void printQudaGaugeParam(QudaGaugeParam *param);
void printQudaInvertParam(QudaInvertParam *param);
void freeGaugeQuda(void);
void loadCloverQuda(void *h_clover, void *h_clovinv,
                      QudaInvertParam *inv_param);
void freeCloverQuda(void);
void invertQuda(void *h_x, void *h_b, QudaInvertParam *param);
void invertMultiShiftQuda(void **_hp_x, void *_hp_b, QudaInvertParam *param);
void dslashQuda(void *h_out, void *h_in, QudaInvertParam *inv_param,
                  QudaParity parity);
void MatQuda(void *h_out, void *h_in, QudaInvertParam *inv_param);
*/

static int isInited = 0;

void
//qhmc_init_quda(int *argc, char ***argv)
qhmc_init_quda(void)
{
  if(isInited) return;
  isInited = 1;

  QudaInitArgs_t init_args;
  //init_args.verbosity = QUDA_SILENT;
  init_args.verbosity = QUDA_SUMMARIZE; /* default */
  //init_args.verbosity = QUDA_VERBOSE;
  //init_args.verbosity = QUDA_DEBUG_VERBOSE;

  const int dim[4] = {nx, ny, nz, nt};
  int status = 0;

  if(is_quda_initialized)return status;

  init_args.layout.device = 0;                                                          // only valid for single-gpu build
  init_args.layout.latsize = dim;
  init_args.layout.machsize = get_logical_dimensions();
  qudaInit(init_args);
}

void
qhmc_fini_quda(void)
{
  endQuda();
}

#if 0
// should take a qopqdp gauge_t and load it
static int
quda_load(lua_State *L)
{
  //void initCommsGridQuda(int nDim, const int *dims, QudaCommsMap func, void *fdata);
  //int device = 0; // FIXME
  //initQuda(device);
  //setVerbosityQuda(QUDA_SUMMARIZE, "QUDA: ", stdout);
  //QudaGaugeParam newQudaGaugeParam(void);
  //void loadGaugeQuda(void *h_gauge, QudaGaugeParam *param);
  return 0;
}
#endif


static int
asqtad_solve_quda(lua_State *L)
{
#define NC QDP_get_nc(qs->field)
  BEGIN_ARGS;
  GET_ASQTAD(a);
  GET_AS_QOPQDP_CVECTOR_ARRAY(nqd, qd);
  GET_QOPQDP_CVECTOR(qs);
  GET_DOUBLE_ARRAY(nm, mass);
  qassert(nqd==nm);
  GET_DOUBLE(resid);

  //int ompnt = omp_get_max_threads();
  //omp_set_num_threads(1);

  QOP_evenodd_t eo=QOP_EVENODD;
  //int prec = 1;
  int max_iter = -1;
  int restart = 500;
  int max_restarts = 5;
  int use_prev_soln = 0;
  double mixed_rsq = 0;

  if(nargs>=nextarg && lua_isstring(L, nextarg)) {
    eo = qopqdp_check_evenodd(L, nextarg);
    nextarg++;
  }

  if(nargs>=nextarg && !lua_isnil(L,nextarg)) {
    if(!lua_istable(L,nextarg)) {
      qlerror0(L,1,"expecting solver paramter table\n");
    }
#define seti(s) lua_getfield(L,nextarg,#s);if(!lua_isnil(L,-1))s=luaL_checkint(L,-1);lua_pop(L,1)
#define setd(s) lua_getfield(L,nextarg,#s);if(!lua_isnil(L,-1))s=luaL_checknumber(L,-1);lua_pop(L,1)
    //seti(prec);
    seti(max_iter);
    seti(restart);
    seti(max_restarts);
    seti(use_prev_soln);
    setd(mixed_rsq);
#undef seti
#undef setd
    nextarg++;
  }
  if(max_iter<0) max_iter = restart*max_restarts;

  END_ARGS;

  //asqtad_set(a, 2);
  //if(h->fatphase || h->longphase) {
  //QOP_asqtad_rephase_field_L_qdp(fla, h->fatphase, h->longphase);
  //}
  QOP_invert_arg_t invarg = QOP_INVERT_ARG_DEFAULT;
  invarg.max_iter = max_iter;
  invarg.restart = restart;
  invarg.max_restarts = max_restarts;
  invarg.evenodd = eo;
  QOP_resid_arg_t resarg[nm], *rap[nm];
  for(int i=0; i<nm; i++) {
    resarg[i] = QOP_RESID_ARG_DEFAULT;
    resarg[i].rsqmin = resid*resid;
    rap[i] = &resarg[i];
  }
  QDP_ColorVector *qqd[nqd];
  for(int i=0; i<nqd; i++) {
    qqd[i] = qd[i]->field;
    if(!use_prev_soln) {
      QDP_V_eq_zero(qqd[i], QDP_all_L(qd[i]->qlat));
    }
  }



  qudaInvert(PRECISION,
             quda_precision, 
             mass,
             inv_args,
             qic->resid,
             qic->relresid,
             fatlink, 
             longlink,
             u0,
             t_src, 
             t_dest,
             &residual,
             &relative_residual, 
             &num_iters);

static struct luaL_Reg quda_reg[] = {
  { "load",           quda_load },
  { NULL, NULL}
};

void
qhmc_open_quda(lua_State *L)
{
  luaL_register(L, "quda", quda_reg);
}
