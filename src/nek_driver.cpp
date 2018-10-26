#include "stream/nek_driver.h"

#include "stream/nek_interface.h"
#include "stream/error.h"

#include <string>

namespace stream {

NekDriver::NekDriver(MPI_Comm comm) : HeatFluidsDriver(comm)
{
  lelg_ = nek_get_lelg();
  lelt_ = nek_get_lelt();
  lx1_ = nek_get_lx1();

  if (active()) {
    MPI_Fint int_comm = MPI_Comm_c2f(comm_.comm);
    C2F_nek_init(static_cast<const int*>(&int_comm));

    nelgt_ = nek_get_nelgt();
    nelt_ = nek_get_nelt();

    init_displs();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

// This version sets the local-to-global element ordering, as ensured by a Gatherv operation.
// It is currently unused, as the coupling does not need to know the local-global ordering.
// void NekDriver::init_mappings() {
//   if(active()) {
//
//     // These arrays are only gathered on the root process
//     if (comm_.rank == 0) {
//       local_counts_.reserve(comm_.size);
//       local_displs_.reserve(comm_.size);
//       local_ordering_.reserve(nelgt_);
//     }
//
//     // Every proc sets a mapping from its local to global element indices.
//     // This mapping is in a Nek5000 common block, but we don't expose it to C++
//     int local_to_global[nelt_];
//     for (int i = 0; i < nelt_; ++i) {
//       local_to_global[i] = nek_get_global_elem(i+1) - 1;
//     }
//
//     // The root proc gets every proc's local element count.
//     comm_.Gather(&nelt_, 1, MPI_INT, local_counts_.data(), 1, MPI_INT);
//
//     if (comm_.rank == 0) {
//       // The root makes a list of data displacements for a Gatherv operation.
//       // Each proc's data will be displaced by the number of local elements on
//       // the previous proc.
//       local_displs_.at(0) = 0;
//       for (int i = 1; i < comm_.rank; ++i) {
//         local_displs_.at(i) = local_displs_.at(i-1) + local_counts_.at(i-1);
//       }
//
//       // The root the gets the local-to-global element mapping for all procs.
//       // This can be used to reorder the data from a Gatherv operation.
//       comm_.Gatherv(
//           local_to_global, nelt_, MPI_INT,
//           local_ordering_.data(), local_counts_.data(), local_displs_.data(), MPI_INT
//       );
//     }
//     else {
//       // Other procs send their local-to-global mapping to the root.
//       comm_.Gatherv(
//           local_to_global, nelt_, MPI_INT,
//           nullptr, nullptr, nullptr, MPI_INT
//       );
//     }
//   }
// }

void NekDriver::init_displs() {
  if(active()) {
    local_counts_.resize(comm_.size);
    local_displs_.resize(comm_.size);

    comm_.Allgather(&nelt_, 1, MPI_INT, local_counts_.data(), 1, MPI_INT);

    local_displs_.at(0) = 0;
    for (int i = 1; i < comm_.size; ++i) {
      local_displs_.at(i) = local_displs_.at(i - 1) + local_counts_.at(i - 1);
    }
  }
}

void NekDriver::init_step() {}

void NekDriver::solve_step() { C2F_nek_solve(); }

void NekDriver::finalize_step() {}

Position NekDriver::get_global_elem_centroid(int global_elem) const
{
  Position centroid;
  err_chk(nek_get_global_elem_centroid(global_elem, &centroid),
      "Could not find centroid of global element " + std::to_string(global_elem));
  return centroid;
}

Position NekDriver::get_local_elem_centroid(int local_elem) const
{
  Position centroid;
  int ierr = nek_get_local_elem_centroid(local_elem, &centroid);
  err_chk(nek_get_local_elem_centroid(local_elem, &centroid),
      "Could not find centroid of local element " + std::to_string(local_elem));
  return centroid;
}

double NekDriver::get_local_elem_volume(int local_elem) const
{
  double volume;
  err_chk(nek_get_local_elem_volume(local_elem, &volume),
      "Could not find volume of local element " + std::to_string(local_elem));
  return volume;
}

double NekDriver::get_local_elem_temperature(int local_elem) const
{
  double temperature;
  err_chk(nek_get_local_elem_temperature(local_elem, &temperature),
      "Could not find temperature of local element " + std::to_string(local_elem));
  return temperature;
}

double NekDriver::get_local_elem_density(int local_elem) const
{
  double density;
  err_chk(nek_get_local_elem_density(local_elem, &density),
          "Could not find density of local element " + std::to_string(local_elem));
  return density;
}

bool NekDriver::global_elem_is_in_rank(int global_elem) const
{
  return (nek_global_elem_is_in_rank(global_elem, comm_.rank) == 1);
}

NekDriver::~NekDriver()
{
  if (active())
    C2F_nek_end();
  MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace stream
