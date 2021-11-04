#include "enrico/surrogate_heat_driver.h"

#include "enrico/vtk_viz.h"
#include "iapws/iapws.h"
#include "openmc/xml_interface.h"
#include "surrogates/heat_xfer_backend.h"
#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xnorm.hpp"
#include "xtensor/xview.hpp"

#include <algorithm> // for fill_n
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iterator> // for back_inserter

namespace enrico {

int ChannelFactory::index_ = 0;
int RodFactory::index_ = 0;

SurrogateHeatDriver::SurrogateHeatDriver(MPI_Comm comm, pugi::xml_node node)
  : HeatFluidsDriver(comm, node)
{
  // Determine thermal-hydraulic parameters for solid phase
  clad_inner_radius_ = node.child("clad_inner_radius").text().as_double();
  clad_outer_radius_ = node.child("clad_outer_radius").text().as_double();
  pellet_radius_ = node.child("pellet_radius").text().as_double();
  n_fuel_rings_ = node.child("fuel_rings").text().as_int();
  n_clad_rings_ = node.child("clad_rings").text().as_int();
  n_pins_x_ = node.child("n_pins_x").text().as_int();
  n_pins_y_ = node.child("n_pins_y").text().as_int();
  n_pins_ = n_pins_x_ * n_pins_y_;
  n_solid_ = n_pins_ * n_axial_ * n_rings() * n_azimuthal_;
  n_fluid_ = n_pins_ * n_axial_;
  pin_pitch_ = node.child("pin_pitch").text().as_double();

  // Determine thermal-hydraulic parameters for fluid phase
  inlet_temperature_ = node.child("inlet_temperature").text().as_double();
  mass_flowrate_ = node.child("mass_flowrate").text().as_double();
  n_channels_ = (n_pins_x_ + 1) * (n_pins_y_ + 1);

  // Determine solver parameters
  if (node.child("max_subchannel_its"))
    max_subchannel_its_ = node.child("max_subchannel_its").text().as_int();
  if (node.child("subchannel_tol_h"))
    subchannel_tol_h_ = node.child("subchannel_tol_h").text().as_double();
  if (node.child("subchannel_tol_p"))
    subchannel_tol_p_ = node.child("subchannel_tol_p").text().as_double();
  if (node.child("heat_tol"))
    heat_tol_ = node.child("heat_tol").text().as_double();

  verbosity_ = verbose::NONE;
  if (node.child("verbosity")) {
    std::string setting = node.child("verbosity").text().as_string();
    if (setting == "none") {
      verbosity_ = verbose::NONE;
    } else if (setting == "low") {
      verbosity_ = verbose::LOW;
    } else if (setting == "high") {
      verbosity_ = verbose::HIGH;
    } else {
      // invalid input for verbosity
      Expects(false);
    }
  }

  // check validity of user input
  Expects(clad_inner_radius_ > 0);
  Expects(clad_outer_radius_ > clad_inner_radius_);
  Expects(pellet_radius_ < clad_inner_radius_);
  Expects(n_fuel_rings_ > 0);
  Expects(n_clad_rings_ > 0);
  Expects(n_pins_x_ > 0);
  Expects(n_pins_y_ > 0);
  Expects(pin_pitch_ > 2.0 * clad_outer_radius_);
  Expects(mass_flowrate_ > 0.0);
  Expects(inlet_temperature_ > 0.0);
  Expects(max_subchannel_its_ > 0);
  Expects(subchannel_tol_h_ > 0.0);
  Expects(subchannel_tol_p_ > 0.0);
  Expects(heat_tol_ > 0.0);

  // Set pin locations, where the center of the assembly is assumed to occur at
  // x = 0, y = 0. It is also assumed that the rod-boundary separation in the
  // x and y directions is the same and equal to half the pitch.
  // TODO: generalize to multi-assembly simulations
  double assembly_width_x = n_pins_x_ * pin_pitch_;
  double assembly_width_y = n_pins_y_ * pin_pitch_;
  double top_left_x = -assembly_width_x / 2.0 + pin_pitch_ / 2.0;
  double top_left_y = assembly_width_y / 2.0 - pin_pitch_ / 2.0;

  pin_centers_.resize({n_pins_, 2});
  for (gsl::index row = 0; row < n_pins_y_; ++row) {
    for (gsl::index col = 0; col < n_pins_x_; ++col) {
      int pin_index = row * n_pins_x_ + col;
      pin_centers_(pin_index, 0) = top_left_x + col * pin_pitch_;
      pin_centers_(pin_index, 1) = top_left_y - row * pin_pitch_;
    }
  }

  // Initialize the channels
  ChannelFactory channel_factory(pin_pitch_, clad_outer_radius_);

  for (std::size_t row = 0; row < n_pins_y_ + 1; ++row) {
    for (std::size_t col = 0; col < n_pins_x_ + 1; ++col) {
      std::size_t a = col / n_pins_x_;
      std::size_t b = row / n_pins_y_;

      if ((row == 0 || row == n_pins_y_) && (col == 0 || col == n_pins_x_))
        channels_.push_back(channel_factory.make_corner(
          {a * (n_pins_x_ - 1) + b * n_pins_x_ * (n_pins_y_ - 1)}));
      else if (row == 0)
        channels_.push_back(channel_factory.make_edge({col - 1, col}));
      else if (row == n_pins_y_)
        channels_.push_back(channel_factory.make_edge(
          {(row - 1) * n_pins_x_ + col - 1, (row - 1) * n_pins_x_ + col}));
      else if (col == 0)
        channels_.push_back(
          channel_factory.make_edge({(row - 1) * n_pins_x_, row * n_pins_x_}));
      else if (col == n_pins_x_)
        channels_.push_back(
          channel_factory.make_edge({row * n_pins_x_ - 1, (row + 1) * n_pins_x_ - 1}));
      else {
        std::size_t i = (row - 1) * n_pins_x_ + col - 1;
        channels_.push_back(
          channel_factory.make_interior({i, i + 1, i + n_pins_x_, i + n_pins_x_ + 1}));
      }
    }
  }

  // Initialize the rods
  RodFactory rod_factory(clad_outer_radius_, clad_inner_radius_, pellet_radius_);
  for (gsl::index rod = 0; rod < n_pins_; ++rod) {
    std::size_t row = rod / n_pins_x_;
    std::size_t col = rod % n_pins_x_;
    std::size_t a = n_pins_x_ + 1;
    rods_.push_back(rod_factory.make_rod(
      {row * a + col, row * a + col + 1, (row + 1) * a + col, (row + 1) * a + col + 1}));
  }

  double total_flow_area = 0.0;
  for (const auto& c : channels_)
    total_flow_area += c.area_;

  channel_flowrates_.resize({n_channels_});
  for (gsl::index i = 0; i < n_channels_; ++i)
    channel_flowrates_(i) = channels_[i].area_ / total_flow_area * mass_flowrate_;

  // Get z values
  z_ = openmc::get_node_xarray<double>(node, "z");
  n_axial_ = z_.size() - 1;

  // Check for visualization input
  if (node.child("viz")) {
    pugi::xml_node viz_node = node.child("viz");
    if (viz_node.attribute("filename")) {
      viz_basename_ = viz_node.attribute("filename").value();
    }

    // if a viz node is found, write final iteration by default
    viz_iterations_ = "final";
    if (viz_node.child("iterations")) {
      viz_iterations_ = viz_node.child("iterations").text().as_string();
    }
    // set other viz values
    if (viz_node.child("resolution")) {
      vtk_radial_res_ = viz_node.child("resolution").text().as_int();
    }
    if (viz_node.child("data")) {
      viz_data_ = viz_node.child("data").text().as_string();
    }
    if (viz_node.child("regions")) {
      viz_regions_ = viz_node.child("regions").text().as_string();
    }
  }

  // Initialize heat transfer solver
  generate_arrays();
};

void SurrogateHeatDriver::generate_arrays()
{
  // Make a radial grid for the clad with equal spacing.
  r_grid_clad_ =
    xt::linspace<double>(clad_inner_radius_, clad_outer_radius_, n_clad_rings_ + 1);

  // Make a radial grid for the fuel with equal spacing.
  r_grid_fuel_ = xt::linspace<double>(0, pellet_radius_, n_fuel_rings_ + 1);

  // Compute the cross-sectional areas of each region
  solid_areas_ = xt::empty<double>({n_rings()});
  for (gsl::index i = 0; i < n_rings(); ++i) {
    if (i < n_fuel_rings_)
      solid_areas_(i) = M_PI * (r_grid_fuel_(i + 1) * r_grid_fuel_(i + 1) -
                                r_grid_fuel_(i) * r_grid_fuel_(i));
    else {
      int r = i - n_fuel_rings_;
      solid_areas_(i) = M_PI * (r_grid_clad_(r + 1) * r_grid_clad_(r + 1) -
                                r_grid_clad_(r) * r_grid_clad_(r));
    }
  }

  if (this->has_coupling_data()) {
    // Create empty arrays for source term and temperature in the solid phase
    source_ = xt::empty<double>({n_pins_, n_axial_, n_rings(), n_azimuthal_});
    solid_temperature_ = xt::empty<double>({n_pins_, n_axial_, n_rings()});

    // Create empty arrays for temperature and density in the fluid phase
    fluid_temperature_ = xt::empty<double>({n_pins_, n_axial_});
    fluid_density_ = xt::empty<double>({n_pins_, n_axial_});
  }
}

int SurrogateHeatDriver::n_local_elem() const
{
  return this->has_coupling_data() ? this->n_global_elem() : 0;
}

std::size_t SurrogateHeatDriver::n_global_elem() const
{
  auto n_solid = n_pins_ * n_axial_ * n_rings() * n_azimuthal_;
  auto n_fluid = n_pins_ * n_axial_;
  return n_solid + n_fluid;
}

std::vector<Position> SurrogateHeatDriver::centroid() const
{
  if (!this->has_coupling_data())
    return {};

  std::vector<Position> centroids;

  // Establish mappings between solid regions and OpenMC cells. The center
  // coordinate for each region in the T/H model is obtained and used to
  // determine the OpenMC cell at that position.
  for (gsl::index i = 0; i < n_pins_; ++i) {
    double x_center = pin_centers_(i, 0);
    double y_center = pin_centers_(i, 1);

    for (gsl::index j = 0; j < n_axial_; ++j) {
      double zavg = 0.5 * (z_(j) + z_(j + 1));

      for (gsl::index k = 0; k < n_rings(); ++k) {
        double ravg;
        if (k < n_fuel_rings_) {
          ravg = 0.5 * (r_grid_fuel_(k) + r_grid_fuel_(k + 1));
        } else {
          int m = k - n_fuel_rings_;
          ravg = 0.5 * (r_grid_clad_(m) + r_grid_clad_(m + 1));
        }

        for (gsl::index m = 0; m < n_azimuthal_; ++m) {
          double m_avg = m + 0.5;
          double theta = 2.0 * m_avg * M_PI / n_azimuthal_;
          double x = x_center + ravg * std::cos(theta);
          double y = y_center + ravg * std::sin(theta);

          // Determine cell instance corresponding to given pin location
          centroids.emplace_back(x, y, zavg);
        }
      }
    }
  }

  // Establish mappings between fluid regions and OpenMC cells; it is assumed that there
  // are no azimuthal divisions in the fluid phase in the OpenMC model so that we
  // can take a point on a 45 degree ray from the pin center. TODO: add a check to make
  // sure that the T/H model is finer than the OpenMC model.

  for (gsl::index i = 0; i < n_pins_; ++i) {
    double x_center = pin_centers_(i, 0);
    double y_center = pin_centers_(i, 1);

    for (gsl::index j = 0; j < n_axial_; ++j) {
      double zavg = 0.5 * (z_(j) + z_(j + 1));
      double l = pin_pitch() / std::sqrt(2.0);
      double d = (l - clad_outer_radius_) / 2.0;
      double x = x_center + (clad_outer_radius_ + d) * std::sqrt(2.0) / 2.0;
      double y = y_center + (clad_outer_radius_ + d) * std::sqrt(2.0) / 2.0;

      // Determine cell instance corresponding to given fluid location
      centroids.emplace_back(x, y, zavg);
    }
  }

  return centroids;
}

std::vector<double> SurrogateHeatDriver::temperature() const
{
  std::vector<double> local_temperatures;

  if (this->has_coupling_data()) {
    for (gsl::index i = 0; i < n_pins_; ++i) {
      for (gsl::index j = 0; j < n_axial_; ++j) {
        for (gsl::index k = 0; k < n_rings(); ++k) {
          for (gsl::index m = 0; m < n_azimuthal_; ++m) {
            local_temperatures.push_back(solid_temperature_(i, j, k));
          }
        }
      }
    }

    for (double T : fluid_temperature_) {
      local_temperatures.push_back(T);
    }
  }

  return local_temperatures;
}

std::vector<double> SurrogateHeatDriver::density() const
{
  std::vector<double> local_densities;

  if (this->has_coupling_data()) {
    // Solid region just gets zeros for densities (not used)
    auto n = n_pins_ * n_axial_ * n_rings() * n_azimuthal_;
    std::fill_n(std::back_inserter(local_densities), n, 0.0);

    // Add fluid densities and return
    for (double rho : fluid_density_) {
      local_densities.push_back(rho);
    }
  }
  return local_densities;
}

int SurrogateHeatDriver::in_fluid_at(int32_t local_elem) const
{
  return local_elem >= n_solid_;
}

std::vector<int> SurrogateHeatDriver::fluid_mask() const
{
  std::vector<int> fluid_mask;

  if (this->has_coupling_data()) {
    auto n_solid = n_pins_ * n_axial_ * n_rings() * n_azimuthal_;
    auto n_fluid = n_pins_ * n_axial_;
    std::fill_n(std::back_inserter(fluid_mask), n_solid, 0);
    std::fill_n(std::back_inserter(fluid_mask), n_fluid, 1);
  }
  return fluid_mask;
}

std::vector<double> SurrogateHeatDriver::volume() const
{
  std::vector<double> volumes;

  if (this->has_coupling_data()) {
    // Volume of solid regions
    for (gsl::index i = 0; i < n_pins_; ++i) {
      for (gsl::index j = 0; j < n_axial_; ++j) {
        double dz = z_(j + 1) - z_(j);
        for (gsl::index k = 0; k < n_rings(); ++k) {
          for (gsl::index m = 0; m < n_azimuthal_; ++m) {
            volumes.push_back(solid_areas_(k) * dz / n_azimuthal_);
          }
        }
      }
    }

    // Volume of fluid regions
    for (gsl::index i = 0; i < n_pins_; ++i) {
      for (gsl::index j = 0; j < n_axial_; ++j) {
        double dz = z_(j + 1) - z_(j);
        double area =
          pin_pitch_ * pin_pitch_ - M_PI * clad_outer_radius_ * clad_outer_radius_;
        volumes.push_back(area * dz);
      }
    }
  }

  return volumes;
}

int SurrogateHeatDriver::set_heat_source_at(int32_t local_elem, double heat)
{
  if (local_elem >= n_pins_ * n_axial_ * n_rings() * n_azimuthal_)
    return 0;

  // Determine indices
  gsl::index pin = local_elem / (n_axial_ * n_rings() * n_azimuthal_);
  gsl::index axial = (local_elem / (n_rings() * n_azimuthal_)) % n_axial_;
  gsl::index ring = (local_elem / n_azimuthal_) % n_rings();
  gsl::index azimuthal = local_elem % n_azimuthal_;

  // Set heat source
  source_(pin, axial, ring, azimuthal) = heat;
  return 0;
}

double SurrogateHeatDriver::rod_axial_node_power(const int pin, const int axial) const
{
  Expects(axial < n_axial_);

  double power = 0.0;
  double dz = z_(axial + 1) - z_(axial);

  for (gsl::index i = 0; i < n_rings(); ++i) {
    for (gsl::index j = 0; j < n_azimuthal_; ++j) {
      power += source_(pin, axial, i, j) * solid_areas_(i) * dz / n_azimuthal_;
    }
  }

  return power;
}

void SurrogateHeatDriver::solve_step()
{
  timer_solve_step.start();
  if (has_coupling_data()) {
    solve_fluid();
    solve_heat();
  }
  timer_solve_step.stop();
}

void SurrogateHeatDriver::solve_fluid()
{
  // determine the power deposition in each channel; the target applications will
  // always be steady-state or pseudo-steady-state cases with no axial conduction such
  // that the power deposition in each channel is independent of a convective heat
  // transfer coefficient and only depends on the rod power at that axial elevation.
  // The channel powers are indexed by channel ID, axial ID
  xt::xtensor<double, 2> channel_powers({n_channels_, n_axial_}, 0.0);
  for (int i = 0; i < n_channels_; ++i) {
    for (int j = 0; j < n_axial_; ++j) {
      for (const auto& rod : channels_[i].rod_ids_)
        channel_powers(i, j) += 0.25 * rod_axial_node_power(rod, j);
    }
  }

  // initial guesses for the fluid solution are uniform temperature (set to the inlet
  // temperature) and uniform pressure  (set to the outlet pressure). These solution
  // fields are defined on channel axial faces. The units used throughout this section
  // are h (kJ/kg), P (MPa), u (m/s), rho (kg/m^3). Unit conversions are performed as
  // necessary on the converged results before being used in the Monte Carlo solver.
  // Enthalpy here requires a factor of 1e-3 to convert from J/kg to kJ/kg.
  xt::xtensor<double, 2> h({n_channels_, n_axial_ + 1},
                           iapws::h1(pressure_bc_, inlet_temperature_));
  xt::xtensor<double, 2> p({n_channels_, n_axial_ + 1}, pressure_bc_);

  // for certain verbosity settings, we will need to save the velocity solutions
  xt::xtensor<double, 2> u = xt::zeros<double>({n_channels_, n_axial_ + 1});

  bool converged = false;
  for (gsl::index iter = 0; iter < max_subchannel_its_; ++iter) {
    // save the previous solution
    xt::xtensor<double, 2> h_old = h;
    xt::xtensor<double, 2> p_old = p;

    // solve each channel independently
    for (gsl::index chan = 0; chan < n_channels_; ++chan) {
      const auto& c = channels_[chan];

      // solve for enthalpy by simple energy balance q = mdot * dh by marching from
      // inlet; divide term on RHS by 1e3 to convert from J/kg to kJ/kg
      h(chan, 0) = iapws::h1(p(chan, 0), inlet_temperature_);
      for (gsl::index axial = 0; axial < n_axial_; ++axial)
        h(chan, axial + 1) =
          h(chan, axial) + 1e-3 * channel_powers(chan, axial) / channel_flowrates_(chan);

      // solve for pressure using one-sided finite difference approximation by
      // marching from outlet and solving the axial momentum equation.
      p(chan, n_axial_) = pressure_bc_;
      for (gsl::index axial = n_axial_; axial > 0; axial--) {
        double rho_high = iapws::rho_from_p_h(p(chan, axial), h(chan, axial));
        double rho_low = iapws::rho_from_p_h(p(chan, axial - 1), h(chan, axial - 1));

        u(chan, axial) = channel_flowrates_(chan) / (rho_high * c.area_);
        u(chan, axial - 1) = channel_flowrates_(chan) / (rho_low * c.area_);

        // factor of 1e-6 needed for convert from Pa to MPa
        p[axial - 1] = p[axial] + 1.0e-6 * (channel_flowrates_(chan) / c.area_ *
                                              (u(chan, axial) - u(chan, axial - 1)) +
                                            g_ * (z_(axial) - z_(axial - 1)) * rho_low);
      }
    }

    // after solving all channels, check for convergence; although all channels are
    // independent, to enable crossflow coupling between channels in the future, this
    // convergence check is performed on all channels together, rather than each
    // separately, since in a more sophisticated solver the channels would all be
    // linked
    auto h_norm = xt::norm_l1(h - h_old)();
    auto p_norm = xt::norm_l1(p - p_old)();

    converged = (h_norm < subchannel_tol_h_) && (p_norm < subchannel_tol_p_);

    if (converged)
      break;

    // check if the solve didn't converge
    if (iter == max_subchannel_its_ - 1) {
      if (verbosity_ >= verbose::LOW) {
        std::cout << "Subchannel solver failed to converge! Enthalpy norm: " << h_norm
                  << " Pressure norm: " << p_norm << std::endl;
      }
    }
  }

  // compute temperature and density from enthalpy and pressure in a cell-centered
  // basis
  xt::xtensor<double, 2> T({n_channels_, n_axial_});
  xt::xtensor<double, 2> rho({n_channels_, n_axial_});

  for (gsl::index chan = 0; chan < n_channels_; ++chan) {
    for (gsl::index axial = 0; axial < n_axial_; ++axial) {
      double h_mean = 0.5 * (h(chan, axial) + h(chan, axial + 1));
      double p_mean = 0.5 * (p(chan, axial) + p(chan, axial + 1));

      T(chan, axial) = iapws::T_from_p_h(p_mean, h_mean);
      rho(chan, axial) = iapws::rho_from_p_h(p_mean, h_mean);
    }
  }

  // After solving the subchannel equations, convert the solution to a rod-centered
  // basis, since this will most likely be the form desired by neutronics codes. At
  // this point only do we apply the conversion of kg/m^3 to g/cm^3 assumed by the
  // neutronics codes.
  for (gsl::index rod = 0; rod < n_pins_; ++rod) {
    for (gsl::index axial = 0; axial < n_axial_; ++axial) {
      fluid_temperature_(rod, axial) = 0.0;
      fluid_density_(rod, axial) = 0.0;

      for (const auto& c : rods_[rod].channel_ids_) {
        fluid_temperature_(rod, axial) += 0.25 * T(c, axial);

        // factor of 1e-3 to convert from kg/m^3 to g/cm^3
        fluid_density_(rod, axial) += 0.25 * rho(c, axial) * 1.0e-3;
      }
    }
  }

  // Perform diagnostic checks if verbosity is sufficiently high
  if (verbosity_ >= verbose::LOW) {
    bool mass_conserved = is_mass_conserved(rho, u);
    bool energy_conserved = is_energy_conserved(rho, u, h, channel_powers);

    Expects(mass_conserved);
    Expects(energy_conserved);
  }
}

bool SurrogateHeatDriver::is_mass_conserved(const xt::xtensor<double, 2>& rho,
                                            const xt::xtensor<double, 2>& u) const
{
  bool mass_conserved = true;

  for (gsl::index axial = 0; axial < n_axial_; ++axial) {
    double mass_flowrate = 0.0;
    for (gsl::index chan = 0; chan < n_channels_; ++chan) {
      double u_cell_centered = 0.5 * (u(chan, axial) + u(chan, axial + 1));
      mass_flowrate += u_cell_centered * channels_[chan].area_ * rho(chan, axial);
    }

    double tol = std::abs(mass_flowrate - mass_flowrate_) / mass_flowrate_;

    if (tol > 1e-3) {
      mass_conserved = false;
    }

    if (verbosity_ == verbose::HIGH) {
      std::cout << "Mass on plane " << axial << " conserved to a tolerance of " << tol
                << std::endl;
    }
  }

  return mass_conserved;
}

bool SurrogateHeatDriver::is_energy_conserved(const xt::xtensor<double, 2>& rho,
                                              const xt::xtensor<double, 2>& u,
                                              const xt::xtensor<double, 2>& h,
                                              const xt::xtensor<double, 2>& q) const
{
  bool energy_conserved = true;

  for (gsl::index axial = 0; axial < n_axial_; ++axial) {
    for (gsl::index chan = 0; chan < n_channels_; ++chan) {
      double u_cell_centered = 0.5 * (u(chan, axial) + u(chan, axial + 1));
      double mass_flowrate = rho(chan, axial) * channels_[chan].area_ * u_cell_centered;

      // conversion factor of 1e3 to convert enthalpy from kJ/kg to J/kg
      double channel_energy_change =
        mass_flowrate * (h(chan, axial + 1) - h(chan, axial)) * 1.0e3;

      double tol = std::abs(channel_energy_change - q(chan, axial)) / q(chan, axial);

      if (tol > 1e-3) {
        energy_conserved = false;
      }

      if (verbosity_ == verbose::HIGH) {
        std::cout << "Energy deposition in channel " << chan << ", axial node " << axial
                  << " conserved to a tolerance of " << tol << std::endl;
      }
    }
  }

  return energy_conserved;
}

void SurrogateHeatDriver::solve_heat()
{
  comm_.message("Solving heat equation...");

  // Convert source to [W/m^3] as expected by Magnolia
  xt::xtensor<double, 3> q = 1e6 * xt::mean(source_, 3);

  // Convert radial grids to [m] as expected by Magnolia
  xt::xtensor<double, 1> r_fuel = 0.01 * r_grid_fuel_;
  xt::xtensor<double, 1> r_clad = 0.01 * r_grid_clad_;

  for (gsl::index i = 0; i < n_pins_; ++i) {
    for (gsl::index j = 0; j < n_axial_; ++j) {
      // approximate cladding surface temperature as equal to the fluid
      // temperature, i.e. this neglects any heat transfer resistance
      double T_co = fluid_temperature_(i, j);

      // Set initial temperature to surface temperature
      xt::view(solid_temperature_, i, j) = T_co;

      solve_steady_nonlin(&q(i, j, 0),
                          T_co,
                          r_fuel.data(),
                          r_clad.data(),
                          n_fuel_rings_,
                          n_clad_rings_,
                          heat_tol_,
                          &solid_temperature_(i, j, 0));
    }
  }
}

double SurrogateHeatDriver::solid_temperature(std::size_t pin,
                                              std::size_t axial,
                                              std::size_t ring) const
{
  return solid_temperature_(pin, axial, ring);
}

double SurrogateHeatDriver::fluid_density(std::size_t pin, std::size_t axial) const
{
  return fluid_density_(pin, axial);
}

double SurrogateHeatDriver::fluid_temperature(std::size_t pin, std::size_t axial) const
{
  return fluid_temperature_(pin, axial);
}

void SurrogateHeatDriver::write_step(int timestep, int iteration)
{
  // TODO: Timers deadlock here...
  // timer_write_step.start();
  if (!has_coupling_data())
    return;

  // if called, but viz isn't requested for the situation,
  // exit early - no output
  if ((iteration < 0 && "final" != viz_iterations_) ||
      (iteration >= 0 && "all" != viz_iterations_)) {
    return;
  }

  // otherwise construct an appropriate filename and write the data
  std::stringstream filename;
  filename << viz_basename_;
  if (iteration >= 0 && timestep >= 0) {
    filename << "_t" << timestep << "_i" << iteration;
  }
  filename << ".vtk";

  SurrogateVtkWriter vtk_writer(*this, vtk_radial_res_, viz_regions_, viz_data_);

  comm_.message("Writing VTK file: " + filename.str());
  vtk_writer.write(filename.str());
  // timer_write_step.stop();
  return;
}

} // namespace enrico
