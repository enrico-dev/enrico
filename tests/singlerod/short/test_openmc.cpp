#include "openmc_driver.h"
#include "mpi.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  {
    stream::OpenmcDriver test_driver(argc, argv, MPI_COMM_WORLD);
    test_driver.init_step();
    test_driver.solve_step();
    test_driver.finalize_step();
  }

  MPI_Finalize();

  return 0;
}
