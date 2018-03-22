#include "drivers.h"
#include "mpi.h"
#include "openmc.h"
#include "nek_interface.h"

// ============================================================================
// OpenMC Driver
// ============================================================================

OpenmcDriver::OpenmcDriver(MPI_Comm comm) : comm(comm) {
  // ROR: 2018-03-22: MPI_Comm_c2f is a macro (in MPICH, at least),
  // so we can't pass something like:
  //     openmc_init(&MPI_Comm_c2f(comm));
  // Hence, the dummy variable.
  MPI_Fint intComm = MPI_Comm_c2f(comm);
  openmc_init(static_cast<const int *>(&intComm));
}

void OpenmcDriver::initStep() {
  openmc_simulation_init();
}

void OpenmcDriver::solveStep() {
  openmc_run();
}

void OpenmcDriver::finalizeStep() {
  openmc_simulation_finalize();
}

OpenmcDriver::~OpenmcDriver() {
  openmc_finalize();
}

// ============================================================================
// Nek5000 Driver
// ============================================================================

NekDriver::NekDriver(MPI_Comm comm) : comm(comm) {
  MPI_Fint intComm = MPI_Comm_c2f(comm);
  C2F_nek_init(static_cast<const int *>(&intComm));
}

void NekDriver::initStep() {
  C2F_nek_init_step();
}

void NekDriver::solveStep() {
  C2F_nek_step();
}

void NekDriver::finalizeStep() {
  C2F_nek_finalize_step();
}

NekDriver::~NekDriver() {
  C2F_nek_end();
}

// ============================================================================
// Coupled Driver
// ============================================================================

CoupledDriver::CoupledDriver(MPI_Comm globalComm, MPI_Comm openmcComm, MPI_Comm nekComm) :
    globalComm(globalComm), openmc(openmcComm), nek(nekComm) {}
