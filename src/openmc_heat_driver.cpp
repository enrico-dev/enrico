#include "stream/openmc_heat_driver.h"

#include "stream/message_passing.h"

#include <gsl/gsl>
#include "openmc/xml_interface.h"
#include "openmc/constants.h"
#include "xtensor/xadapt.hpp"
#include "xtensor/xio.hpp"

#include <cmath>
#include <unordered_map>

namespace stream {

SurrogateHeatDriver::SurrogateHeatDriver(MPI_Comm comm, pugi::xml_node node)
  : HeatFluidsDriver{comm}
{
  // Determine heat transfer solver parameters
  double clad_ir = node.child("clad_inner_radius").text().as_double();
  double clad_or = node.child("clad_outer_radius").text().as_double();
  double pellet = node.child("pellet_radius").text().as_double();
  int fuel_rings = node.child("fuel_rings").text().as_int();
  int clad_rings = node.child("clad_rings").text().as_int();

  // Get pin locations
  // TODO: Switch to get_node_xarray on OpenMC update
  auto pin_locations = openmc::get_node_array<double>(node, "pin_centers");
  if (pin_locations.size() % 2 != 0) {
    throw std::runtime_error{"Length of <pin_centers> must be a multiple of two"};
  }

  // Convert to xtensor
  auto npins = pin_locations.size() / 2;
  std::vector<std::size_t> shape {npins, 2};
  pin_centers_ = xt::adapt(pin_locations, shape);

  // Get z values
  // TODO: Switch to get_node_xarray on OpenMC update
  auto z_values = openmc::get_node_array<double>(node, "z");
  z_ = xt::adapt(z_values);

  // Initialize heat transfer solver
  solver_ = std::make_unique<HeatSolver>(clad_ir, clad_or, pellet, npins);
  solver_->set_rings(fuel_rings, clad_rings);
};

void SurrogateHeatDriver::set_heat_source()
{
};

OpenmcHeatDriver::OpenmcHeatDriver(MPI_Comm comm, pugi::xml_node node)
  : comm_{comm}
{
  // Get internode and cross-node communicators
  MPI_Comm intranode;
  MPI_Comm crossnode;
  get_node_comms(comm_.comm, 1, &intranode, &crossnode);

  // Get power in [W]
  power_ = node.child("power").text().as_double();

  // Initialize OpenMC and surrogate heat drivers
  openmc_driver_ = std::make_unique<OpenmcDriver>(crossnode);
  heat_driver_ = std::make_unique<SurrogateHeatDriver>(comm, node);

  // Save internode communicator
  intranode_comm_ = Comm{intranode};

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
      int nrings = heat_driver_->solver_->n_fuel_rings_;
      const auto& radii = heat_driver_->solver_->r_grid_fuel_;
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
      mats.push_back(c.material_index_);
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
  // Determine heat source based on OpenMC tally results
  auto Q = openmc_driver_->heat_source(power_);

  auto nz = heat_driver_->z_.size() - 1;
  auto npins = heat_driver_->pin_centers_.shape()[0];
  int ring_index = 0;
  int pin_index = 0;
  for (int i = 0; i < npins; ++i) {
    for (int j = 0; j < nz; ++ j) {
      // Loop over radial rings
      int nrings = heat_driver_->solver_->n_fuel_rings_;
      for (int k = 0; k < nrings; ++k) {
        // Determine cell instances present in this ring
        const auto& cell_instances = ring_to_cell_inst_[ring_index];
        auto n = cell_instances.size();

        // Get average Q value across each azimuthal segment
        double q_avg = 0.0;
        for (int k = 0; k < n; ++k) {
          q_avg += Q(cell_instances[k]);
        }
        q_avg /= n;

        // Set Q in appropriate (pin, ring)
        heat_driver_->solver_->source_(pin_index, k) = q_avg;
        ++ring_index;
      }
      ++pin_index;
    }
  }
}

void OpenmcHeatDriver::update_temperature()
{

}

} // namespace stream
