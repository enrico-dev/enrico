#include "openmc/capi.h"
#include "openmc_driver.h"
#include "error.h"

namespace stream {

OpenmcDriver::OpenmcDriver(int argc, char* argv[], MPI_Comm comm)
    : TransportDriver(comm)
{
  if (active()) {
    err_chk(openmc_init(argc, argv, &comm), openmc_err_msg);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void OpenmcDriver::init_step() { err_chk(openmc_simulation_init(), openmc_err_msg); }

void OpenmcDriver::solve_step() { err_chk(openmc_run(), openmc_err_msg); }

void OpenmcDriver::finalize_step() { err_chk(openmc_simulation_finalize(), openmc_err_msg); }

OpenmcDriver::~OpenmcDriver()
{
  if (active())
    err_chk(openmc_finalize(), openmc_err_msg);
  MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace stream
