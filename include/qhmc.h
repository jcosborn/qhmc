#include "qhmc_base.h"

void qhmc_init_libs(int *argc, char ***argv);
void qhmc_fini_libs(void);
void qhmc_open_libs(lua_State* L);

int luaopen_bc(lua_State *L);

int luaopen_lfs (lua_State *L);

#ifdef QHMC_HAVE_QOPQDP
#include "qhmc_qopqdp.h"
#define qhmc_master() qhmc_qopqdp_master()
#else
#define qhmc_master() 1
#endif

void qhmc_init_quda(int *argc, char ***argv);
void qhmc_fini_quda(void);
void qhmc_open_quda(lua_State* L);
