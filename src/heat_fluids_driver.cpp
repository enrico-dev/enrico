#include "enrico/heat_fluids_driver.h"
#include "gsl/gsl"

namespace enrico {

HeatFluidsDriver::HeatFluidsDriver(MPI_Comm comm, double pressure_bc)
  : Driver(comm)
  , pressure_bc_(pressure_bc)
{
  Expects(pressure_bc_ > 0.0);
}

void HeatFluidsDriver::init_displs()
{
  if (active()) {
    local_counts_.resize(comm_.size);
    local_displs_.resize(comm_.size);

    int32_t n_local = this->n_local_elem();
    comm_.Allgather(&n_local, 1, MPI_INT32_T, local_counts_.data(), 1, MPI_INT32_T);

    local_displs_.at(0) = 0;
    for (gsl::index i = 1; i < comm_.size; ++i) {
      local_displs_.at(i) = local_displs_.at(i - 1) + local_counts_.at(i - 1);
    }
  }
}

std::vector<Position> HeatFluidsDriver::centroids() const
{
  // Get local centroids on each rank
  auto local_centroids = this->centroid_local();

  // Gather local centroids onto root process
  return this->gather(local_centroids);
}

std::vector<double> HeatFluidsDriver::volumes() const
{
  // Get local volumes on each rank
  auto local_volumes = this->volume_local();

  // Gather local centroids onto root process
  return this->gather(local_volumes);
}

}
