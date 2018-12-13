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

using namespace xt::placeholders;

namespace stream {

class SurrogateHeatDriver : public HeatFluidsDriver {

enum class VTKData {heat = 0, neutronics, all};

public:
  //! Initializes heat-fluids surrogate with the given MPI communicator.
  //!
  //! \param comm  The MPI communicator used to initialze the surrogate
  //! \param node  XML node containing settings for surrogate
  explicit SurrogateHeatDriver(MPI_Comm comm, pugi::xml_node node);

  //! Finalizes heat-fluids surrogate.
  ~SurrogateHeatDriver() { };

  //! Initializes timestep for heat-fluids surrogate solver
  void init_step() { };

  //! Solves the heat equation at a single timestep in each specified region
  void solve_step();

  //! Finalizes the timestep for the heat-fluids surrogate solver
  void finalize_step() { };

  //! Number of rings in fuel and clad
  //! \return Number of rings
  std::size_t n_rings() { return n_fuel_rings_ + n_clad_rings_; }

  //! Write data to VTK
  void to_vtk(std::string filename = "magnolia.vtk");

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

private:
  //! Create internal arrays used for heat equation solver
  void generate_arrays();

  class VisualizationPin {
    // some constant values
    const int WEDGE_TYPE_ = 13;
    const int WEDGE_SIZE_ = 6;
    const int HEX_TYPE_ = 12;
    const int HEX_SIZE_ = 8;
    const int INVALID_CONN_ = -1;

  public:
    VisualizationPin(double x, double y,
                     xt::xtensor<double, 1> z_grid,
                     xt::xtensor<double, 1> r_grid,
                     xt::xtensor<double, 1> c_grid,
                     int t_res);

    // methods
    xt::xtensor<double, 2> create_ring(double radius, int t_resolution);
    xt::xtensor<double, 3> pin_points();
    xt::xtensor<int, 4> pin_connectivity();

    xt::xtensor<double, 3> clad_points();
    xt::xtensor<int, 4> clad_connectivity();

    inline int num_points() { return points_per_plane_ * (axial_divs_ + 1); }
    inline int num_entries() { return (WEDGE_SIZE_ + 1) * (axial_divs_ * t_res_) +
        (HEX_SIZE_ + 1) * t_res_ * (radial_divs_ - 1) * axial_divs_; }
    inline int num_cells() { return axial_divs_ * radial_divs_ * t_res_; }
    inline int conn_entry_size() { return HEX_SIZE_ + 1; }

    // members
    double x_, y_;
    int t_res_;
    int cells_per_plane_;
    int points_per_plane_;
    int radial_divs_;
    int axial_divs_;
    xt::xtensor<double, 1> z_grid_;
    xt::xtensor<double, 1> r_grid_;
    xt::xtensor<double, 1> c_grid_;
    xt::xtensor<double, 1> data_;
  };
}; // end SurrogateHeatDriver

} // namespace stream

#endif // STREAM_HEAT_DRIVER_H
