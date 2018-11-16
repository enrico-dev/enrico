#include "stream/openmc_heat_driver.h"

#include "stream/message_passing.h"

#include <gsl/gsl>
#include "openmc/constants.h"

#include <cmath>
#include <unordered_map>

namespace stream {

OpenmcHeatDriver::OpenmcHeatDriver(MPI_Comm comm, pugi::xml_node node)
  : comm_{comm}
{
  // Get power in [W]
  power_ = node.child("power").text().as_double();

  // Initialize OpenMC and surrogate heat drivers
  openmc_driver_ = std::make_unique<OpenmcDriver>(comm);
  pugi::xml_node surr_node = node.child("heat_surrogate");
  heat_driver_ = std::make_unique<SurrogateHeatDriver>(comm, surr_node);

  // Create mappings for fuel pins and setup tallies for OpenMC
  init_mappings();
  init_tallies();

}

void OpenmcHeatDriver::init_mappings()
{
  using openmc::PI;

  const auto& centers = heat_driver_->pin_centers_;
  const auto& z = heat_driver_->z_;
  auto nz = z.size() - 1;
  auto npins = centers.shape()[0];

  std::unordered_map<int32_t, int> tracked;
  int ring_index = 0;
  // TODO: Don't hardcode number of azimuthal segments
  int n_azimuthal = 4;
  for (int i = 0; i < npins; ++i) {
    for (int j = 0; j < nz; ++j) {
      // Get average z value
      double zavg = 0.5*(z(j) + z(j + 1));

      // Loop over radial rings
      int nrings = heat_driver_->n_fuel_rings_;
      const auto& radii = heat_driver_->r_grid_fuel_;
      for (int k = 0; k < nrings; ++k) {
        double ravg = 0.5*(radii(k) + radii(k + 1));

        for (int m = 0; m < n_azimuthal; ++m) {
          double theta = 2.0*m*PI/n_azimuthal + 0.01;
          double x = ravg * std::cos(theta);
          double y = ravg * std::sin(theta);

          // Determine cell instance corresponding to given pin location
          Position r {x, y, zavg};
          CellInstance c {r};
          if (tracked.find(c.material_index_) == tracked.end()) {
            openmc_driver_->cells_.push_back(c);
            tracked[c.material_index_] = openmc_driver_->cells_.size() - 1;
          }

          // Map OpenMC material to ring and vice versa
          int32_t array_index = tracked[c.material_index_];
          cell_inst_to_ring_[array_index].push_back(ring_index);
          ring_to_cell_inst_[ring_index].push_back(array_index);
        }

      ++ring_index;
      }
    }
  }
}

void OpenmcHeatDriver::init_tallies()
{
  if (openmc_driver_->active()) {
    // Build vector of material indices
    std::vector<int32_t> mats;
    for (const auto& c : openmc_driver_->cells_) {
      mats.push_back(c.material_index_ - 1);
    }
    openmc_driver_->create_tallies(mats);
  }
}

void OpenmcHeatDriver::solve_step()
{
  // Solve neutron transport
  if (openmc_driver_->active()) {
    openmc_driver_->solve_step();
  }
  comm_.Barrier();

  // Update heat source
  update_heat_source();

  // Solve heat equation
  if (heat_driver_->active()) {
    heat_driver_->solve_step();
  }
  comm_.Barrier();

  // Update temperature in OpenMC
  update_temperature();
}

void OpenmcHeatDriver::update_heat_source()
{
  // zero out heat source
  for (auto& val : heat_driver_->source_) {
    val = 0.0;
  }

  // Determine heat source based on OpenMC tally results
  auto Q = openmc_driver_->heat_source(power_);

  int ring_index = 0;
  for (int i = 0; i < heat_driver_->n_pins_; ++i) {
    for (int j = 0; j < heat_driver_->n_axial_; ++j) {
      // Loop over radial rings
      int nrings = heat_driver_->n_fuel_rings_;
      for (int k = 0; k < nrings; ++k) {
        // Get average Q value across each azimuthal segment
        const auto& cell_instances = ring_to_cell_inst_[ring_index];
        double q_avg = 0.0;
        for (auto idx : cell_instances) {
          q_avg += Q(idx);
        }
        q_avg /= cell_instances.size();

        // Set Q in appropriate (pin, axial, ring)
        heat_driver_->source_.at(i, j, k) = q_avg;
        ++ring_index;
      }
    }
  }
}

void OpenmcHeatDriver::update_temperature()
{

}

} // namespace stream
