//! \file heat_driver.h
//! Driver for Magnolia's heat transfer solver
#ifndef ENRICO_SURROGATE_HEAT_DRIVER_H
#define ENRICO_SURROGATE_HEAT_DRIVER_H

#include "enrico/geom.h"
#include "enrico/heat_fluids_driver.h"

#include <gsl/gsl>
#include <mpi.h>
#include <pugixml.hpp>
#include <xtensor/xtensor.hpp>

#include <cstddef>

namespace enrico {

//! Struct containing geometric information for a flow channel
struct Channel {
  //! Channel index
  int index_;

  //! Channel flow area
  double area_;

  //! Vector of rod IDs connected to this channel, all with a fractional perimeter
  //! in contact with the channel equal to 0.25
  std::vector<std::size_t> rod_ids_;
};

//! Struct containing geometry information for a cylindrical solid rod
struct Rod {
  //! Rod index
  int index_;

  //! Rod cladding outer radius
  double clad_outer_radius_;

  //! Rod cladding inner radius
  double clad_inner_radius_;

  //! Rod pellet radius
  double pellet_radius_;

  //! Vector of channel IDs connected to this rod, all with a fractional perimeter
  //! in contact with the rod equal to 0.25
  std::vector<std::size_t> channel_ids_;
};

//! Class to construct flow channels for a Cartesian lattice of pins
class ChannelFactory {
public:
  ChannelFactory(double pitch, double rod_radius)
    : pitch_(pitch)
    , radius_(rod_radius)
    , interior_area_(pitch_ * pitch_ - M_PI * radius_ * radius_)
  {}

  //! Make a corner subchannel connected to given rods
  Channel make_corner(const std::vector<std::size_t>& rods) const
  {
    Channel c;
    c.index_ = index_++;
    c.area_ = 0.25 * interior_area_;
    c.rod_ids_ = rods;
    return c;
  }

  //! Make an edge subchannel connected to given rods
  Channel make_edge(const std::vector<std::size_t>& rods) const
  {
    Channel c;
    c.index_ = index_++;
    c.area_ = 0.5 * interior_area_;
    c.rod_ids_ = rods;
    return c;
  }

  //! Make an interior subchannel connected to given rods
  Channel make_interior(const std::vector<std::size_t>& rods) const
  {
    Channel c;
    c.index_ = index_++;
    c.area_ = interior_area_;
    c.rod_ids_ = rods;
    return c;
  }

private:
  //! rod pitch
  double pitch_;

  //! rod outer radius
  double radius_;

  //! interior channel flow area, which is proportional to flow areas for all
  //! other channel types
  double interior_area_;

  //! index of constructed channel
  static int index_;
};

class RodFactory {
public:
  RodFactory(double clad_OR, double clad_IR, double pellet_OR)
    : clad_outer_r_(clad_OR)
    , clad_inner_r_(clad_IR)
    , pellet_outer_r_(pellet_OR)
  {}

  //! Make a rod connected to given channels
  Rod make_rod(const std::vector<std::size_t>& channels) const
  {
    Rod r;
    r.index_ = index_++;
    r.clad_inner_radius_ = clad_inner_r_;
    r.clad_outer_radius_ = clad_outer_r_;
    r.pellet_radius_ = pellet_outer_r_;
    r.channel_ids_ = channels;
    return r;
  }

private:
  //! Cladding outer radius
  double clad_outer_r_;

  //! Cladding inner radius
  double clad_inner_r_;

  //! Pellet outer radius
  double pellet_outer_r_;

  //! Index of constructed rod
  static int index_;
};

/**
 * Class providing surrogate thermal-hydraulic solution for a Cartesian
 * bundle of rods with upwards-flowing coolant. A conduction model is used
 * for the solid phase, with axial conduction neglected. The solid phase is
 * linked to the fluid phase by conjugate heat transfer, which is treated
 * here with a pseudo-steady-state approach where the power entering the fluid
 * matches the power in the rod at that axial elevation. It is assumed that
 * there is zero thermal resistance between the rod and the fluid.
 *
 * The fluid solution is obtained with a very simplified "subchannel" method
 * that neglects all crossflow terms between channels such that the method
 * is more akin to a single-pin analysis in a coolant-centered basis. The
 * enthalpy is solved by simply axial energy balance, while the axial momentum
 * equation is solved for pressure (the mass flow rate in each channel being
 * fixed) while neglecting friction effects.
 */
class SurrogateHeatDriver : public HeatFluidsDriver {
public:
  //! Initializes heat-fluids surrogate with the given MPI communicator.
  //!
  //! \param comm  The MPI communicator used to initialze the surrogate
  //! \param node  XML node containing settings for surrogate
  SurrogateHeatDriver(MPI_Comm comm, pugi::xml_node node);

  //! Verbosity options for printing simulation results
  enum class verbose { NONE, LOW, HIGH };

  bool has_coupling_data() const final { return comm_.rank == 0; }

  //! Get the number of local mesh elements
  //! \return Number of local mesh elements
  int n_local_elem() const override;

  //! Get the number of global mesh elements
  //! \return Number of global mesh elements
  std::size_t n_global_elem() const override;

  int in_fluid_at(int32_t local_elem) const override;

  //! Set the heat source for a given local element
  //!
  //! \param local_elem A local element ID
  //! \param heat A heat source term
  //! \return Error code
  int set_heat_source_at(int32_t local_elem, double heat) override;

  //! Solves the heat-fluids surrogate solver
  void solve_step() final;

  void solve_heat();

  void solve_fluid();

  //! Returns Number of rings in fuel and clad
  std::size_t n_rings() const { return n_fuel_rings_ + n_clad_rings_; }

  //! Returns cladding inner radius
  double clad_inner_radius() const { return clad_inner_radius_; }

  //! Returns cladding outer radius
  double clad_outer_radius() const { return clad_outer_radius_; }

  //! Returns pellet outer radius
  double pellet_radius() const { return pellet_radius_; }

  //! Returns number of fuel rings
  std::size_t n_fuel_rings() const { return n_fuel_rings_; }

  //! Returns number of clad rings
  std::size_t n_clad_rings() const { return n_clad_rings_; }

  //! Returns number of pins in x-direction
  std::size_t n_pins_x() const { return n_pins_x_; }

  //! Returns number of pins in y-direction
  std::size_t n_pins_y() const { return n_pins_y_; }

  //! Returns number of solid elements
  std::size_t n_solid_;

  //! Returns number of fluid elements
  std::size_t n_fluid_;

  //! Returns pin pitch
  double pin_pitch() const { return pin_pitch_; }

  //! Returns inlet temperature boundary condition in [K]
  double inlet_temperature() const { return inlet_temperature_; }

  //! Returns inlet mass flowrate boundary condition in [kg/s]
  double mass_flowrate() const { return mass_flowrate_; }

  //! Returns maximum number of subchannel iterations
  std::size_t max_subchannel_its() const { return max_subchannel_its_; }

  //! Returns subchannel convergence tolerance for enthalpy
  double subchannel_tol_h() const { return subchannel_tol_h_; }

  //! Returns subchannel convergence tolerance for pressure
  double subchannel_tol_p() const { return subchannel_tol_p_; }

  //! Returns convergence tolerance for solid energy equation
  double heat_tol() const { return heat_tol_; }

  //! Write data to VTK
  void write_step(int timestep, int iteration) final;

  //! Returns solid temperature in [K] for given region
  double solid_temperature(std::size_t pin, std::size_t axial, std::size_t ring) const;

  //! Returns fluid density in [g/cm^3] for given region
  double fluid_density(std::size_t pin, std::size_t axial) const;

  //! Returns fluid temperature in [K] for given region
  double fluid_temperature(std::size_t pin, std::size_t axial) const;

  // Data on fuel pins
  xt::xtensor<double, 2> pin_centers_; //!< (x,y) values for center of fuel pins
  xt::xtensor<double, 1> z_;           //!< Bounding z-values for axial segments
  std::size_t n_axial_;                //!< number of axial segments
  std::size_t n_azimuthal_{4};         //!< number of azimuthal segments

  //! Total number of pins
  std::size_t n_pins_;

  // Dimensions for a single fuel pin axial segment
  double clad_outer_radius_;     //!< clad outer radius in [cm]
  double clad_inner_radius_;     //!< clad inner radius in [cm]
  double pellet_radius_;         //!< fuel pellet radius in [cm]
  std::size_t n_fuel_rings_{20}; //!< number of fuel rings
  std::size_t n_clad_rings_{2};  //!< number of clad rings

  //!< Channels in the domain
  std::vector<Channel> channels_;

  //!< Rods in the domain
  std::vector<Rod> rods_;

  //! Mass flowrate for coolant-centered channels; this is determine by distributing
  //! a total inlet mass flowrate among the channels based on the fractional flow area.
  xt::xtensor<double, 1> channel_flowrates_;

  // solver variables and settings
  xt::xtensor<double, 4>
    source_; //!< heat source for each (pin, axial segment, ring, azimuthal segment)
  xt::xtensor<double, 1> r_grid_clad_; //!< radii of each clad ring in [cm]
  xt::xtensor<double, 1> r_grid_fuel_; //!< radii of each fuel ring in [cm]

  //! Cross-sectional areas of rings in fuel and cladding
  xt::xtensor<double, 1> solid_areas_;

  // visualization
  std::string viz_basename_{
    "heat_surrogate"}; //!< base filename for visualization files (default: magnolia)
  std::string viz_iterations_{
    "none"};                    //!< visualization iterations to write (none, all, final)
  std::string viz_data_{"all"}; //!< visualization data to write
  std::string viz_regions_{"all"}; //!< visualization regions to write
  size_t vtk_radial_res_{20};      //!< radial resolution of resulting vtk files

private:
  //! Get temperature of local mesh elements
  //! \return Temperature of local mesh elements in [K]
  std::vector<double> temperature() const override;

  //! Get density of local mesh elements
  //! \return Density of local mesh elements in [g/cm^3]
  std::vector<double> density() const override;

  //! States whether each local region is in fluid
  //! \return For each local region, 1 if region is in fluid and 0 otherwise
  std::vector<int> fluid_mask() const override;

  //! Get centroids of local mesh elements
  //! \return Centroids of local mesh elements
  std::vector<Position> centroid() const override;

  //! Get volumes of local mesh elements
  //! \return Volumes of local mesh elements
  std::vector<double> volume() const override;

  //! Create internal arrays used for heat equation solver
  void generate_arrays();

  //! Channel index in terms of row, column index
  int channel_index(int row, int col) const { return row * (n_pins_x_ + 1) + col; }

  //! Rod power at a given node in a given pin, computed by integrating the heat source
  //! (assumed constant in each ring) over the pin.
  //! \param pin   pin index
  //! \param axial axial index
  double rod_axial_node_power(const int pin, const int axial) const;

  //! Diagnostic function to assess whether the mass is conserved by the subchannel
  //! solver by comparing the mass flowrate in each axial plane (at cell-centered
  //! positions) to the specified inlet mass flowrate.
  //! \param rho density in a cell-centered basis
  //! \param u   axial velocity in a face-centered basis
  bool is_mass_conserved(const xt::xtensor<double, 2>& rho,
                         const xt::xtensor<double, 2>& u) const;

  //! Diagnostic function to assess whether the energy is conserved by the subchannel
  //! solver by comparing the energy deposition in each channel in each axial plane
  //! (at cell-centered positions) to the powers of the rods connected to that channel.
  //! \param rho density in a cell-centered basis
  //! \param u   axial velocity in a face-centered basis
  //! \param h   enthalpy in a face-centered basis
  //! \param q   powers in each channel in a cell-centered basis
  bool is_energy_conserved(const xt::xtensor<double, 2>& rho,
                           const xt::xtensor<double, 2>& u,
                           const xt::xtensor<double, 2>& h,
                           const xt::xtensor<double, 2>& q) const;

  //!< solid temperature in [K] for each (pin, axial segment, ring)
  xt::xtensor<double, 3> solid_temperature_;

  //! Flow areas for coolant-centered channels
  xt::xtensor<double, 1> channel_areas_;

  //! Fluid temperature in a rod-centered basis indexed by rod ID and axial ID
  xt::xtensor<double, 2> fluid_temperature_;

  //! Fluid density in [g/cm^3] in a rod-centered basis indexed by rod ID and axial ID
  xt::xtensor<double, 2> fluid_density_;

  //! Number of pins in the x-direction in a Cartesian grid
  std::size_t n_pins_x_;

  //! Number of pins in the y-direction in a Cartesian grid
  std::size_t n_pins_y_;

  //! Pin pitch, assumed the same for the x and y directions
  double pin_pitch_;

  //! Inlet fluid temperature [K]
  double inlet_temperature_;

  //! Mass flowrate of fluid into the domain [kg/s]
  double mass_flowrate_;

  //! Number of channels
  std::size_t n_channels_;

  //! Maximum number of iterations for subchannel solution, set to a default value
  //! of 100 if not set by the user
  int max_subchannel_its_ = 100;

  //! Convergence tolerance on enthalpy for the subchannel solution for use in
  //! convergence based on the L-1 norm, set to a default value of 1e-2
  double subchannel_tol_h_ = 1e-2;

  //! Convergence tolerance on pressure for the subchannel solution for use in
  //! convergence based on the L-1 norm, set to a default value of 1e-2
  double subchannel_tol_p_ = 1e-2;

  //! Convergence tolerance for solid temperature solution, set to a default value
  //! of 1e-4
  double heat_tol_ = 1e-4;

  //! Gravitational acceleration
  const double g_ = 9.81;

  //! Verbosity setting for printing simulation results; defaults to NONE
  verbose verbosity_ = verbose::NONE;

}; // end SurrogateHeatDriver

} // namespace enrico

#endif // ENRICO_SURROGATE_HEAT_DRIVER_H
