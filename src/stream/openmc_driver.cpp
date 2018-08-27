#include "openmc/capi.h"
#include "openmc_driver.h"

namespace stream {

OpenmcDriver::OpenmcDriver(int argc, char *argv[], MPI_Comm comm)
    : TransportDriver(comm) {
  if (active()) {
    openmc_init(argc, argv, &comm);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void OpenmcDriver::init_step() { openmc_simulation_init(); }

void OpenmcDriver::solve_step() { openmc_run(); }

void OpenmcDriver::finalize_step() { openmc_simulation_finalize(); }

OpenmcDriver::~OpenmcDriver() {
  if (active())
    openmc_finalize();
  MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace stream
