#include "stream/base_drivers.h"
#include "stream/message_passing.h"
#include "stream/openmc_nek_driver.h"
#include <mpi.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  {
    stream::OpenmcNekDriver test_driver(1.0, MPI_COMM_WORLD);

    if (test_driver.openmc_driver_->active()) {
      test_driver.openmc_driver_->init_step();
      test_driver.openmc_driver_->solve_step();
      test_driver.openmc_driver_->finalize_step();
    }
    test_driver.comm_.Barrier();

    test_driver.update_heat_source();

    if (test_driver.nek_driver_->active()) {
      test_driver.nek_driver_->init_step();
      test_driver.nek_driver_->solve_step();
      test_driver.nek_driver_->finalize_step();
    }
    test_driver.comm_.Barrier();

    test_driver.update_temperature();
  }

  MPI_Finalize();
}
