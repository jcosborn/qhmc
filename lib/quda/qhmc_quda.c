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

void
qhmc_init_quda(int *argc, char ***argv)
{
  int device = 0; // FIXME
  initQuda(device);
  setVerbosityQuda(QUDA_SUMMARIZE, "QUDA: ", stdout);
}

void
qhmc_fini_quda(void)
{
  endQuda();
}

static int
quda_load(lua_State *L)
{
  //QudaGaugeParam newQudaGaugeParam(void);
  //void loadGaugeQuda(void *h_gauge, QudaGaugeParam *param);
  return 0;
}

static struct luaL_Reg quda_reg[] = {
  { "load",           quda_load },
  { NULL, NULL}
};

void
qhmc_open_quda(lua_State *L)
{
  luaL_register(L, "quda", quda_reg);
}
