#include "enrico/openmc_driver.h"
#include <mpi.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  {
    enrico::OpenmcDriver test_driver(MPI_COMM_WORLD);
    test_driver.init_step();
    test_driver.solve_step(0);
    test_driver.finalize_step();
  }

  MPI_Finalize();

  return 0;
}
