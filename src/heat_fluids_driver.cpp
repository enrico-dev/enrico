#include "enrico/heat_fluids_driver.h"

#include <gsl/gsl>
#include <openmc/vendor/pugixml/src/pugixml.hpp>
#include <xtensor/xadapt.hpp>

namespace enrico {

HeatFluidsDriver::HeatFluidsDriver(MPI_Comm comm, pugi::xml_node node)
  : Driver(comm)
{
  pressure_bc_ = node.child("pressure_bc").text().as_double();
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

xt::xtensor<double, 1> HeatFluidsDriver::temperature() const
{
  // Get local tempratures on each rank
  auto local_temperatures = this->temperature_local();

  // Gather all the local element temperatures onto the root
  auto global_temperatures = this->gather(local_temperatures);

  // only the return value from root should be used, or else a broadcast added here
  return xt::adapt(global_temperatures);
}

xt::xtensor<double, 1> HeatFluidsDriver::density() const
{
  // Get local densities on each rank
  auto local_densities = this->density_local();

  // Gather all local element densities onto the root
  auto global_densities = this->gather(local_densities);

  return xt::adapt(global_densities);
}

std::vector<int> HeatFluidsDriver::fluid_mask() const
{
  // Get local fluid masks
  auto local_fluid_mask = this->fluid_mask_local();

  // Gather all the local fluid masks onto the root
  return this->gather(local_fluid_mask);
}

}
