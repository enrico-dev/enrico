#include "drivers.h"
#include "message_passing.h"
#include "mpi.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // openmcComm is split from MPI_COMM_WORLD.  It will contain 1 proc per node.
  MPI_Comm openmcComm = MPI_COMM_NULL;
  getInternodeSubComm(MPI_COMM_WORLD, 1, &openmcComm);

  MPI_Comm coupledComm = MPI_COMM_WORLD;
  MPI_Comm nekComm = MPI_COMM_WORLD;

  auto *testDriver =
      new OpenmcNekDriver(argc, argv, coupledComm, openmcComm, nekComm);

  testDriver->openmcDriver.initStep();
  testDriver->openmcDriver.solveStep();
  testDriver->openmcDriver.finalizeStep();

  testDriver->nekDriver.initStep();
  testDriver->nekDriver.solveStep();
  testDriver->nekDriver.finalizeStep();

  delete testDriver;

  MPI_Finalize();
}