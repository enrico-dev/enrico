#include "drivers.h"
#include "mpi.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  {
    stream::OpenmcDriver testDriver(argc, argv, MPI_COMM_WORLD);
    testDriver.init_step();
    testDriver.solve_step();
    testDriver.finalize_step();
  }

  MPI_Finalize();

  return 0;
}
