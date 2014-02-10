#include "qhmc_config.h"
#include "qhmc.h"
#include "qhmc_main.h"

int
main(int argc, char *argv[])
{
  qhmc_init_libs(&argc, &argv);

  lmain(argc, argv);

  qhmc_fini_libs();

  return 0;
}
