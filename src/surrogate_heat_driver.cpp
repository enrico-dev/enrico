#include "enrico/surrogate_heat_driver.h"

#include "enrico/vtk_viz.h"
#include "openmc/xml_interface.h"
#include "surrogates/heat_xfer_backend.h"
#include "iapws/iapws.h"
#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xnorm.hpp"

#include <iostream>

namespace enrico {

int ChannelFactory::index_ = 0;

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
    subchannel_tol_h_= node.child("subchannel_tol_h").text().as_double();
  if (node.child("subchannel_tol_p"))
    subchannel_tol_p_= node.child("subchannel_tol_p").text().as_double();
  if (node.child("heat_tol"))
    heat_tol_ = node.child("heat_tol").text().as_double();

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
  double assembly_width_x = n_pins_x_ * pin_pitch_;
  double assembly_width_y = n_pins_y_ * pin_pitch_;

  pin_centers_.resize({n_pins_});
  for (int row = 0; row < n_pins_y_; ++row) {
    for (int col = 0; col < n_pins_x_; ++col) {
      int pin_index = row * n_pins_x_ + col;
      pin_centers_(pin_index, 0) = -assembly_width_x / 2.0 + pin_pitch_ / 2.0 + col * pin_pitch_;
      pin_centers_(pin_index, 1) = assembly_width_y / 2.0 - (pin_pitch_ / 2.0 + row * pin_pitch_);
    }
  }

  // Initialize the channels
  ChannelFactory channel_factory(pin_pitch_, clad_outer_radius_);

  for (int row = 0; row < n_pins_y_ + 1; ++row) {
    for (int col = 0; col < n_pins_x_ + 1; ++col) {
      int a = col / n_pins_x_;
      int b = row / n_pins_y_;

      if ((row == 0 || row == n_pins_y_) && (col == 0 || col == n_pins_x_))
        channels_.push_back(channel_factory.make_corner({a * (n_pins_x_ - 1) + b * n_pins_x_ * (n_pins_y_ - 1)}));
      else if (row == 0)
        channels_.push_back(channel_factory.make_edge({col - 1, col}));
      else if (row == n_pins_y_)
        channels_.push_back(channel_factory.make_edge({(row - 1) * n_pins_x_ + col - 1, (row - 1) * n_pins_x_ + col}));
      else if (col == 0)
        channels_.push_back(channel_factory.make_edge({(row - 1) * n_pins_x_, row * n_pins_x_}));
      else if (col == n_pins_x_)
        channels_.push_back(channel_factory.make_edge({row * n_pins_x_ - 1, (row + 1) * n_pins_x_ - 1}));
      else {
        int i = (row - 1) * n_pins_x_ + col - 1;
        channels_.push_back(channel_factory.make_interior({i, i + 1, i + n_pins_x_, i + n_pins_x_ + 1}));
      }
    }
  }

  double total_flow_area = 0.0;
  for (const auto c : channels_)
    total_flow_area += c.area_;

  channel_flowrates_.resize({n_channels_});
  for (int i = 0; i < n_channels_; ++i)
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
  for (int i = 0; i < n_rings(); ++i) {
    if (i < n_fuel_rings_)
      solid_areas_(i) = M_PI * (r_grid_fuel_(i + 1) * r_grid_fuel_(i + 1) - r_grid_fuel_(i) * r_grid_fuel_(i));
    else {
      int r = i - n_fuel_rings_;
      solid_areas_(i) = M_PI * (r_grid_clad_(r + 1) * r_grid_clad_(r + 1) - r_grid_clad_(r) * r_grid_clad_(r));
    }
  }

  // Create empty arrays for source term and temperature
  source_ = xt::empty<double>({n_pins_, n_axial_, n_rings()});
  temperature_ = xt::empty<double>({n_pins_, n_axial_, n_rings()});
  density_ = xt::zeros<double>({n_pins_, n_axial_, n_rings()});
  fluid_mask_ = xt::zeros<int>({n_pins_, n_axial_, n_rings()});
}

void SurrogateHeatDriver::solve_step()
{
  solve_heat();
  solve_fluid();
}

void SurrogateHeatDriver::solve_fluid()
{
  // convenience function for determining the rod power in a given axial layer
  auto rod_axial_node_power = [this](int pin, int axial) {
    Expects(axial <= n_axial_);

    double power = 0.0;
    double dz = z_(axial + 1) - z_(axial);

    for (int i = 0; i < n_rings(); ++i)
      power += source_(pin, axial, i) * solid_areas_(i) * dz;

    return power;
  };

  // determine the power deposition in each channel; the target applications will always
  // be steady-state or pseudo-steady-state cases with no axial conduction such that
  // the power deposition in each channel is independent of a convective heat transfer
  // coefficient and only depends on the rod power at that axial elevation. The channel
  // powers are indexed by channel ID, axial ID
  xt::xtensor<double, 2> channel_powers = xt::zeros<double>({n_channels_, n_axial_});
  for (int i = 0; i < n_channels_; ++i) {
    for (int j = 0; j < n_axial_; ++j) {
      for (const auto& rod : channels_[i].rod_ids_)
        channel_powers(i, j) += 0.25 * rod_axial_node_power(rod, j);
    }
  }

  // initial guesses for the fluid solution are uniform temperature (set to the inlet
  // temperature) and uniform pressure  (set to the outlet pressure). These solution
  // fields are defined on channel axial faces
  xt::xtensor<double, 2> h = xt::zeros<double>({n_channels_, n_axial_ + 1});
  xt::xtensor<double, 2> p = xt::zeros<double>({n_channels_, n_axial_ + 1});
  std::fill(h.begin(), h.end(), iapws::h1(pressure_bc_, inlet_temperature_));
  std::fill(p.begin(), p.end(), pressure_bc_);

  bool converged = false;
  for (int iter = 0; iter < max_subchannel_its_; ++iter) {
    // save the previous solution
    xt::xtensor<double, 2> h_old = h;
    xt::xtensor<double, 2> p_old = p;

    // solve each channel independently
    for (int chan = 0; chan < n_channels_; ++chan) {
      const auto& c = channels_[chan];

      // solve for enthalpy by simple energy balance q = mdot * dh by marching from inlet
      h(chan, 0) = iapws::h1(p(chan, 0), inlet_temperature_);
      for (int axial = 0; axial < n_axial_; ++axial)
        h(chan, axial + 1) = h(chan, axial) + channel_powers(chan, axial) / channel_flowrates_(chan);

      // solve for pressure using one-sided finite difference approximation by marching
      // from outlet and solving the axial momentum equation.
      p(chan, n_axial_) = pressure_bc_;
      for (int axial = n_axial_; axial > 0; axial--) {
        double rho_high = 1.0e-3 / iapws::nu1(p[axial], iapws::T_from_p_h(p[axial], h[axial]));
        double rho_low = 1.0e-3 / iapws::nu1(p[axial - 1], iapws::T_from_p_h(p[axial - 1], h[axial - 1]));

        double u_high = channel_flowrates_(chan) / (rho_high * c.area_);
        double u_low = channel_flowrates_(chan) / (rho_low * c.area_);

        p[axial - 1] = p[axial] + channel_flowrates_(chan) / c.area_ * (u_high - u_low) +
          g_ * (z_(axial) - z_(axial - 1)) * rho_low;
      }
    }

    // after solving all channels, check for convergence; although all channels are independent,
    // to enable crossflow coupling between channels in the future, this convergence check is
    // performed on all channels together, rather than each separately, since in a more sophisticated
    // solver the channels would all be linked
    auto n_expr_h = xt::norm_l1(h - h_old);
    double h_norm = n_expr_h();
    auto n_expr_p = xt::norm_l1(p - p_old);
    double p_norm = n_expr_p();

    converged = (h_norm < subchannel_tol_h_) && (p_norm < subchannel_tol_p_);

    if (converged)
      break;
  }
}

void SurrogateHeatDriver::solve_heat()
{
  std::cout << "Solving heat equation...\n";
  // NuScale inlet temperature
  double T_co = 523.15;

  // Set initial temperature
  for (auto& T : temperature_) {
    T = T_co;
  }

  // Convert source to [W/m^3] as expected by Magnolia
  xt::xtensor<double, 3> q = 1e6 * source_;

  // Convert radial grids to [m] as expected by Magnolia
  xt::xtensor<double, 1> r_fuel = 0.01 * r_grid_fuel_;
  xt::xtensor<double, 1> r_clad = 0.01 * r_grid_clad_;

  for (int i = 0; i < n_pins_; ++i) {
    for (int j = 0; j < n_axial_; ++j) {
      solve_steady_nonlin(&q(i, j, 0),
                          T_co,
                          r_fuel.data(),
                          r_clad.data(),
                          n_fuel_rings_,
                          n_clad_rings_,
                          heat_tol_,
                          &temperature_(i, j, 0));
    }
  }
}

xt::xtensor<double, 1> SurrogateHeatDriver::temperature() const
{
  return xt::flatten(temperature_);
}

double SurrogateHeatDriver::temperature(int pin, int axial, int ring) const
{
  return temperature_(pin, axial, ring);
}

xt::xtensor<double, 1> SurrogateHeatDriver::density() const
{
  return xt::flatten(density_);
}

xt::xtensor<int, 1> SurrogateHeatDriver::fluid_mask() const
{
  return xt::flatten(fluid_mask_);
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

  std::cout << "Writing VTK file: " << filename.str() << "\n";
  vtk_writer.write(filename.str());
  return;
}

} // namespace enrico
