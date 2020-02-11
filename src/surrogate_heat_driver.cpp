#include "enrico/surrogate_heat_driver.h"

#include "enrico/vtk_viz.h"
#include "iapws/iapws.h"
#include "openmc/xml_interface.h"
#include "surrogates/heat_xfer_backend.h"
#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xnorm.hpp"
#include "xtensor/xview.hpp"

#include <iostream>

namespace enrico {

int ChannelFactory::index_ = 0;
int RodFactory::index_ = 0;

SurrogateHeatDriver::SurrogateHeatDriver(MPI_Comm comm,
                                         double pressure_bc,
                                         pugi::xml_node node)
  : HeatFluidsDriver(comm, pressure_bc)
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

  for (gsl::index row = 0; row < n_pins_y_ + 1; ++row) {
    for (gsl::index col = 0; col < n_pins_x_ + 1; ++col) {
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
  // TODO: Switch to get_node_xarray on OpenMC update
  auto z_values = openmc::get_node_array<double>(node, "z");
  z_ = xt::adapt(z_values);
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

  // Create empty arrays for source term and temperature in the solid phase
  source_ = xt::empty<double>({n_pins_, n_axial_, n_rings()});
  solid_temperature_ = xt::empty<double>({n_pins_, n_axial_, n_rings()});

  // Create empty arrays for temperature and density in the fluid phase
  fluid_temperature_ = xt::empty<double>({n_pins_, n_axial_});
  fluid_density_ = xt::empty<double>({n_pins_, n_axial_});
}

double SurrogateHeatDriver::rod_axial_node_power(const int pin, const int axial) const
{
  Expects(axial < n_axial_);

  double power = 0.0;
  double dz = z_(axial + 1) - z_(axial);

  for (gsl::index i = 0; i < n_rings(); ++i)
    power += source_(pin, axial, i) * solid_areas_(i) * dz;

  return power;
}

void SurrogateHeatDriver::solve_step()
{
  solve_fluid();
  solve_heat();
}

void SurrogateHeatDriver::solve_fluid()
{
  // determine the power deposition in each channel; the target applications will always
  // be steady-state or pseudo-steady-state cases with no axial conduction such that
  // the power deposition in each channel is independent of a convective heat transfer
  // coefficient and only depends on the rod power at that axial elevation. The channel
  // powers are indexed by channel ID, axial ID
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

      // solve for enthalpy by simple energy balance q = mdot * dh by marching from inlet;
      // divide term on RHS by 1e3 to convert from J/kg to kJ/kg
      h(chan, 0) = iapws::h1(p(chan, 0), inlet_temperature_);
      for (gsl::index axial = 0; axial < n_axial_; ++axial)
        h(chan, axial + 1) =
          h(chan, axial) + 1e-3 * channel_powers(chan, axial) / channel_flowrates_(chan);

      // solve for pressure using one-sided finite difference approximation by marching
      // from outlet and solving the axial momentum equation.
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
    // separately, since in a more sophisticated solver the channels would all be linked
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

  // compute temperature and density from enthalpy and pressure in a cell-centered basis
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

  // After solving the subchannel equations, convert the solution to a rod-centered basis,
  // since this will most likely be the form desired by neutronics codes. At this point
  // only do we apply the conversion of kg/m^3 to g/cm^3 assumed by the neutronics codes.
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
  xt::xtensor<double, 3> q = 1e6 * source_;

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

xt::xtensor<double, 1> SurrogateHeatDriver::temperature() const
{
  xt::xarray<double> Ts = {xt::flatten(solid_temperature_)};
  xt::xarray<double> Tf = {xt::flatten(fluid_temperature_)};
  return xt::concatenate(xt::xtuple(Ts, Tf), 0);
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

xt::xtensor<double, 1> SurrogateHeatDriver::density() const
{
  return xt::flatten(fluid_density_);
}

void SurrogateHeatDriver::write_step(int timestep, int iteration)
{
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
  return;
}

} // namespace enrico
