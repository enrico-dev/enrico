#include "enrico/openmc_driver.h"
#include <mpi.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  {
    enrico::OpenmcDriver test_driver(enrico::Comm(MPI_COMM_WORLD));
    test_driver.init_step();
    test_driver.solve_step();
    test_driver.write_step(0, 0);
    test_driver.finalize_step();
  }

  MPI_Finalize();

  return 0;
}
