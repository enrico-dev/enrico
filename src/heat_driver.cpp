#include "stream/heat_driver.h"

#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xstrided_view.hpp"
#include "heat_xfer_backend.h"
#include "openmc/xml_interface.h"

namespace stream {

SurrogateHeatDriver::SurrogateHeatDriver(MPI_Comm comm, pugi::xml_node node)
  : HeatFluidsDriver{comm}
{
  // Determine heat transfer solver parameters
  clad_inner_radius_ = node.child("clad_inner_radius").text().as_double();
  clad_outer_radius_ = node.child("clad_outer_radius").text().as_double();
  pellet_radius_ = node.child("pellet_radius").text().as_double();
  n_fuel_rings_ = node.child("fuel_rings").text().as_int();
  n_clad_rings_ = node.child("clad_rings").text().as_int();

  // Get pin locations
  // TODO: Switch to get_node_xarray on OpenMC update
  auto pin_locations = openmc::get_node_array<double>(node, "pin_centers");
  if (pin_locations.size() % 2 != 0) {
    throw std::runtime_error{"Length of <pin_centers> must be a multiple of two"};
  }

  // Convert to xtensor
  n_pins_ = pin_locations.size() / 2;
  std::vector<std::size_t> shape {n_pins_, 2};
  pin_centers_ = xt::adapt(pin_locations, shape);

  // Get z values
  // TODO: Switch to get_node_xarray on OpenMC update
  auto z_values = openmc::get_node_array<double>(node, "z");
  z_ = xt::adapt(z_values);
  n_axial_ = z_.size() - 1;

  // Initialize heat transfer solver
  generate_arrays();
};

void SurrogateHeatDriver::generate_arrays()
{
  // Make a radial grid for the fuel with equal spacing.
  r_grid_clad_ = xt::linspace<double>(
    clad_inner_radius_,
    clad_outer_radius_,
    n_clad_rings_ + 1
  );

  // Make a radial grid for the clad with equal spacing.
  r_grid_fuel_ = xt::linspace<double>(
    0, pellet_radius_, n_fuel_rings_ + 1);

  // Create empty arrays for source term and temperature
  auto total_axial = n_pins_ * n_axial_;
  source_ = xt::empty<double>({total_axial, n_rings()});
  temperature_ = xt::empty<double>({total_axial, n_rings()});
}

void SurrogateHeatDriver::set_heat_source(gsl::span<double> q)
{
  auto source_flat = xt::flatten(source_);
  for (int i = 0; i < q.size(); ++i) {
    source_flat(i) = q[i];
  }
}

void SurrogateHeatDriver::solve_step(gsl::span<double> T_co)
{
  for (int i = 0; i < n_pins_; ++i) {
    solve_steady_nonlin(&source_(i), T_co[i], r_grid_fuel_.data(), r_grid_clad_.data(),
      n_fuel_rings_, n_clad_rings_, tol_, &temperature_(i));
  }
}

} // namespace stream
