//! \file heat_driver.h
//! Driver for Magnolia's heat transfer solver
#ifndef ENRICO_HEAT_DRIVER_H
#define ENRICO_HEAT_DRIVER_H

#include "driver.h"

#include <gsl/gsl>
#include <mpi.h>
#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

#include <cstddef>

namespace enrico {

class SurrogateHeatDriver : public Driver {
public:
  //! Initializes heat-fluids surrogate with the given MPI communicator.
  //!
  //! \param comm  The MPI communicator used to initialze the surrogate
  //! \param node  XML node containing settings for surrogate
  explicit SurrogateHeatDriver(MPI_Comm comm, pugi::xml_node node);

  //! Solves the heat-fluids surrogate solver
  virtual void solve_step() override;

  //! Number of rings in fuel and clad
  //! \return Number of rings
  std::size_t n_rings() { return n_fuel_rings_ + n_clad_rings_; }

  //! Write data to VTK
  virtual void write_step(int timestep = -1, int iteration = -1) override;

  // Data on fuel pins
  xt::xtensor<double, 2> pin_centers_; //!< (x,y) values for center of fuel pins
  xt::xtensor<double, 1> z_;           //!< Bounding z-values for axial segments
  std::size_t n_pins_;                 //!< number of fuel pins
  std::size_t n_axial_;                //!< number of axial segments

  // Dimensions for a single fuel pin axial segment
  double clad_outer_radius_;      //!< clad outer radius in [cm]
  double clad_inner_radius_;      //!< clad inner radius in [cm]
  double pellet_radius_;          //!< fuel pellet radius in [cm]
  std::size_t n_fuel_rings_ {20}; //!< number of fuel rings
  std::size_t n_clad_rings_ {2};  //!< number of clad rings

  // solver variables and settings
  double tol_;                         //!< tolerance on convergence
  xt::xtensor<double, 3> source_;      //!< heat source for each (axial segment, ring)
  xt::xtensor<double, 3> temperature_; //!< temperature in [K] for each (axial segment, ring)
  xt::xtensor<double, 1> r_grid_clad_; //!< radii of each clad ring in [cm]
  xt::xtensor<double, 1> r_grid_fuel_; //!< radii of each fuel ring in [cm]

  // visualization
  std::string viz_basename_{"heat_surrogate"}; //!< base filename for visualization files (default: magnolia)
  std::string viz_iterations_{"none"};   //!< visualization iterations to write (none, all, final)
  std::string viz_data_{"all"};          //!< visualization data to write
  std::string viz_regions_{"all"};       //!< visualization regions to write
  size_t vtk_radial_res_{20};            //!< radial resolution of resulting vtk files

private:
  //! Create internal arrays used for heat equation solver
  void generate_arrays();
}; // end SurrogateHeatDriver

} // namespace enrico

#endif // ENRICO_HEAT_DRIVER_H
