#include "drivers.h"
#include "message_passing.h"
#include "mpi.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // openmc_comm is split from MPI_COMM_WORLD.  It will contain 1 proc per node.
  MPI_Comm openmc_comm = MPI_COMM_NULL;
  stream::get_internode_sub_comm(MPI_COMM_WORLD, 1, &openmc_comm);

  MPI_Comm coupled_comm = MPI_COMM_WORLD;
  MPI_Comm nek_comm = MPI_COMM_WORLD;

  {
    stream::OpenmcNekDriver test_driver(argc, argv, coupled_comm, openmc_comm, nek_comm);

    if (test_driver.openmc_driver_.proc_info_.comm != MPI_COMM_NULL) {
      test_driver.openmc_driver_.init_step();
      test_driver.openmc_driver_.solve_step();
      test_driver.openmc_driver_.finalize_step();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (test_driver.nek_driver_.proc_info_.comm != MPI_COMM_NULL) {
      test_driver.nek_driver_.init_step();
      test_driver.nek_driver_.solve_step();
      test_driver.nek_driver_.finalize_step();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Finalize();
}