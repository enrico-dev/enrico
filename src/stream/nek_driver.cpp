#include "nek_driver.h"

namespace stream {

NekDriver::NekDriver(MPI_Comm comm) : HeatFluidsDriver(comm) {
  lelg_ = nek_get_lelg();
  lelt_ = nek_get_lelt();
  lx1_ = nek_get_lx1();

  if (active()) {
    MPI_Fint int_comm = MPI_Comm_c2f(comm_.comm);
    C2F_nek_init(static_cast<const int *>(&int_comm));

    nelgt_ = nek_get_nelgt();
    nelt_ = nek_get_nelt();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void NekDriver::init_step() {}

void NekDriver::solve_step() { C2F_nek_solve(); }

void NekDriver::finalize_step() {}

Position NekDriver::get_global_elem_centroid(int global_elem) const {
  Position centroid;
  int ierr = nek_get_global_elem_centroid(global_elem, &centroid);
  return centroid;
}

NekDriver::~NekDriver() {
  if (active())
    C2F_nek_end();
  MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace stream
