#include "qhmc_common.h"

#ifdef HAVE_LFS
#include "src/lfs.h"
#else
#define luaopen_lfs(L)
#endif

#ifdef HAVE_QOPQDP
#include "../qopqdp/qhmc_qopqdp.h"
#endif
