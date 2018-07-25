#include "drivers.h"
#include "message_passing.h"
#include "mpi.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // openmcComm is split from MPI_COMM_WORLD.  It will contain 1 proc per node.
  MPI_Comm openmcComm = MPI_COMM_NULL;
  stream::get_internode_sub_comm(MPI_COMM_WORLD, 1, &openmcComm);

  MPI_Comm coupledComm = MPI_COMM_WORLD;
  MPI_Comm nekComm = MPI_COMM_WORLD;

  {
    stream::OpenmcNekDriver testDriver(argc, argv, coupledComm, openmcComm, nekComm);

    if (testDriver.openmc_driver_.proc_info_.comm != MPI_COMM_NULL) {
      testDriver.openmc_driver_.init_step();
      testDriver.openmc_driver_.solve_step();
      testDriver.openmc_driver_.finalize_step();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (testDriver.nek_driver_.proc_info_.comm != MPI_COMM_NULL) {
      testDriver.nek_driver_.init_step();
      testDriver.nek_driver_.solve_step();
      testDriver.nek_driver_.finalize_step();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Finalize();
}