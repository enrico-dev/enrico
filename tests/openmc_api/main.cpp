#include "drivers.h"
#include "mpi.h"

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  auto * testDriver = new OpenmcDriver(MPI_COMM_WORLD);
  testDriver->initStep();
  testDriver->solveStep();
  testDriver->finalizeStep();
  delete testDriver;

  // HDF5 cleanup occurs in OpenmcDriver's destructor, so MPI_Finalize may have already been called
  int isFinalized;
  MPI_Finalized(&isFinalized);
  if (!isFinalized)
    MPI_Finalize();

  return 0;
}
