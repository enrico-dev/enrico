#include "drivers.h"
#include "mpi.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  {
    NekDriver testDriver(MPI_COMM_WORLD);
    testDriver.initStep();
    testDriver.solveStep();
    testDriver.finalizeStep();
  }

  MPI_Finalize();

  return 0;
}
