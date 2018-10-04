//! \file heat_driver.h
//! Driver for Magnolia's heat transfer solver
#ifndef STREAM_HEAT_DRIVER_H
#define STREAM_HEAT_DRIVER_H

#include "base_drivers.h"

#include <gsl/gsl>
#include <mpi.h>
#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

#include <cstddef>

namespace stream {

class SurrogateHeatDriver : public HeatFluidsDriver {
public:
  explicit SurrogateHeatDriver(MPI_Comm comm, pugi::xml_node node);
  ~SurrogateHeatDriver() { };

  void init_step() { };
  void solve_step();
  void finalize_step() { };

  void set_heat_source(gsl::span<double> q);
  std::size_t n_rings() { return n_fuel_rings_ + n_clad_rings_; }

  // Data on fuel pins
  xt::xtensor<double, 2> pin_centers_; //!< (x,y) values for center of fuel pins
  xt::xtensor<double, 1> z_;  //!< Bounding z-values for axial segments
  std::size_t n_pins_;                //!< number of fuel pins
  std::size_t n_axial_;               //!< number of axial segments

  // Dimensions for a single fuel pin axial segment
  double clad_outer_radius_;  //!< clad outer radius in [m]
  double clad_inner_radius_;  //!< clad inner radius in [m]
  double pellet_radius_;      //!< fuel pellet radius in [m]
  std::size_t n_fuel_rings_ {20};     //!< number of fuel rings
  std::size_t n_clad_rings_ {2};      //!< number of clad rings

  // solver variables and settings
  double tol_; //!< tolerance on convergence
  xt::xtensor<double, 3> source_;      //!< heat source for each (axial segment, ring)
  xt::xtensor<double, 3> temperature_; //!< temperature in [K] for each (axial segment, ring)
  xt::xtensor<double, 1> r_grid_clad_; //!< radii of each clad ring in [cm]
  xt::xtensor<double, 1> r_grid_fuel_; //!< radii of each fuel ring in [cm]

private:
  void generate_arrays();
};

} // namespace stream

#endif // STREAM_HEAT_DRIVER_H
