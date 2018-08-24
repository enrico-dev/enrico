#include "drivers.h"
#include "mpi.h"
#include "nek_interface.h"
#include "openmc/capi.h"
#include "stream_geom.h"

#include <algorithm> // for max, fill, copy
#include <iterator> // for back_inserter
#include <unordered_set>

namespace stream {

// ============================================================================
// HeatFluids Driver
// ============================================================================

bool HeatFluidsDriver::active() const {
  return comm_.comm != MPI_COMM_NULL;
}

// ============================================================================
// Transport Driver
// ============================================================================

bool TransportDriver::active() const {
  return comm_.comm != MPI_COMM_NULL;
}

// ============================================================================
// OpenMC Driver
// ============================================================================

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
