#include "drivers.h"
#include "mpi.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  {
    stream::NekDriver test_driver(MPI_COMM_WORLD);
    test_driver.init_step();
    test_driver.solve_step();
    test_driver.finalize_step();
  }

  MPI_Finalize();

  return 0;
}
