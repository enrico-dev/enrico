//! \file heat_driver.h
//! Interface to Magnolia's heat transfer solver
#ifndef STREAM_HEAT_DRIVER_H
#define STREAM_HEAT_DRIVER_H

#include <cstddef> // for size_t

#include <gsl/gsl>
#include "xtensor/xtensor.hpp"

namespace stream {

class HeatSolver {
public:
  HeatSolver(double clad_outer, double clad_inner, double pellet,
    int n_pins);

  void solve_step(gsl::span<double> T_co);

  void set_heat_source(gsl::span<double> q);

  void set_radii(double clad_outer, double clad_inner, double pellet);
  void set_rings(int n_fuel, int n_clad);
  void set_num_pins(int n);
  int n_rings() { return n_fuel_rings_ + n_clad_rings_; }

  // fuel pin dimensions
  double clad_outer_radius_;   //!< clad outer radius in [m]
  double clad_inner_radius_;   //!< clad inner radius in [m]
  double pellet_radius_;       //!< fuel pellet radius in [m]
  int n_fuel_rings_ {20}; //!< number of fuel rings
  int n_clad_rings_ {2}; //!< number of clad rings
  int n_pins_; //!< number of fuel pins

  // solver variables and settings
  double tol_; //!< tolerance on convergence
  xt::xtensor<double, 2> source_; //!< heat source for each pin / ring
  xt::xtensor<double, 2> temperature_; //!< temperature in [K] for each pin / ring
  xt::xtensor<double, 1> r_grid_clad_;
  xt::xtensor<double, 1> r_grid_fuel_;
private:
  void generate_arrays();
};

} // namespace stream

#endif // STREAM_HEAT_DRIVER_H
