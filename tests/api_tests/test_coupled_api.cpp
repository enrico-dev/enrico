#include "drivers.h"
#include "message_passing.h"
#include "mpi.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // openmcComm is split from MPI_COMM_WORLD.  It will contain 1 proc per node.
  MPI_Comm openmcComm = MPI_COMM_NULL;
  stream::getInternodeSubComm(MPI_COMM_WORLD, 1, &openmcComm);

  MPI_Comm coupledComm = MPI_COMM_WORLD;
  MPI_Comm nekComm = MPI_COMM_WORLD;

  {
    stream::OpenmcNekDriver testDriver(argc, argv, coupledComm, openmcComm, nekComm);

    if (testDriver.openmcDriver.procInfo.comm != MPI_COMM_NULL) {
      testDriver.openmcDriver.initStep();
      testDriver.openmcDriver.solveStep();
      testDriver.openmcDriver.finalizeStep();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (testDriver.nekDriver.procInfo.comm != MPI_COMM_NULL) {
      testDriver.nekDriver.initStep();
      testDriver.nekDriver.solveStep();
      testDriver.nekDriver.finalizeStep();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Finalize();
}