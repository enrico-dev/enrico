#include "drivers.h"
#include "mpi.h"

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  auto *testDriver = new OpenmcDriver(MPI_COMM_WORLD);
  testDriver->initStep();
  testDriver->solveStep();
  testDriver->finalizeStep();
  delete testDriver;

  MPI_Finalize();

  return 0;
}
