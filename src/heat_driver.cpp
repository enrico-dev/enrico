#include "stream/heat_driver.h"

#include "xtensor/xbuilder.hpp"
#include "xtensor/xstrided_view.hpp"
#include "heat_xfer_backend.h"

namespace stream {

HeatSolver::HeatSolver(double clad_outer, double clad_inner, double pellet, int n_pins)
  : n_pins_{n_pins}, clad_outer_radius_{clad_outer},
    clad_inner_radius_{clad_inner}, pellet_radius_{pellet}
{
  generate_arrays();
}

void HeatSolver::set_radii(double clad_outer, double clad_inner, double pellet)
{
  clad_outer_radius_ = clad_outer;
  clad_inner_radius_ = clad_inner;
  pellet_radius_ = pellet;
  generate_arrays();
}

void HeatSolver::set_rings(int n_fuel, int n_clad)
{
  n_fuel_rings_ = n_fuel;
  n_clad_rings_ = n_clad;
  generate_arrays();
}

void HeatSolver::set_num_pins(int n)
{
  n_pins_ = n;
  generate_arrays();
}

void HeatSolver::generate_arrays()
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
  source_ = xt::empty<double>({n_pins_, n_rings()});
  temperature_ = xt::empty<double>({n_pins_, n_rings()});
}

void HeatSolver::set_heat_source(gsl::span<double> q)
{
  auto source_flat = xt::flatten(source_);
  for (int i = 0; i < q.size(); ++i) {
    source_flat(i) = q[i];
  }
}

void HeatSolver::solve_step(gsl::span<double> T_co)
{
  for (int i = 0; i < n_pins_; ++i) {
    solve_steady_nonlin(&source_(i), T_co[i], r_grid_fuel_.data(), r_grid_clad_.data(),
      n_fuel_rings_, n_clad_rings_, tol_, &temperature_(i));
  }
}

} // namespace stream
