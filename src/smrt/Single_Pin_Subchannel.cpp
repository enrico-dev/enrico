#include <cmath>
#include <iostream>

#include "smrt/Single_Pin_Subchannel.h"

// vendored includes
#include "openmc/string_utils.h"
#include <gsl/gsl>
#include <iapws.h>

namespace enrico {
//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
Single_Pin_Subchannel::Single_Pin_Subchannel(RCP_PL& parameters,
                                             const std::vector<double>& delta_z)
  : d_delta_z(delta_z)
  , d_T_inlet(-1.0)
  , d_p_exit(-1.0)
{
  // Read parameters from parameter list
  d_tol = parameters->get("tolerance", 1.0e-6);
  d_max_iters = parameters->get("max_iters", 100);

  // Convert dz from cm to m
  for (auto& val : d_delta_z)
    val *= 1e-2;

  auto verb = parameters->get("verbosity", std::string("none"));
  openmc::to_lower(verb);
  if (verb == "none") {
    d_verbosity = NONE;
  } else if (verb == "low") {
    d_verbosity = LOW;
  } else if (verb == "high") {
    d_verbosity = HIGH;
  } else {
    // Unrecognized verbosity.
    Expects(false);
  }
}

//---------------------------------------------------------------------------//
// Solve for temperature and density in subchannel
//
// Note that the subchannel solver internally works with density in units
// of kg/m^3 (to facilitate use of water material property correlations),
// but is converted to g/cm^3 upon output.
//---------------------------------------------------------------------------//
void Single_Pin_Subchannel::solve(const std::vector<double>& power,
                                  std::vector<double>& temperature,
                                  std::vector<double>& density)
{
  Expects(power.size() == d_delta_z.size());
  Expects(temperature.size() == d_delta_z.size());
  Expects(density.size() == d_delta_z.size());
  Expects(d_T_inlet > 0.0);
  Expects(d_p_exit > 0.0);

  int num_regions = d_delta_z.size();

  // Initial guesses for h and p on edges of axial regions. Note that enthalpy
  // is in units of kJ/kg and pressure is in units of MPa as needed for IAPWS
  // functions
  std::vector<double> h(num_regions + 1, iapws::h1(d_p_exit, d_T_inlet));
  std::vector<double> p(num_regions + 1, d_p_exit);

  // Gravitational constant
  constexpr double g = 9.81;

  // Convert W to kW
  constexpr double w_to_kw = 1e-3;

  bool converged = false;
  int it;
  double h_diff, p_diff;
  for (it = 1; it <= d_max_iters; ++it) {
    std::vector<double> h_old = h;
    std::vector<double> p_old = p;

    // Solve for enthalpy from inlet to outlet
    h[0] = iapws::h1(p[0], d_T_inlet);
    for (int k = 0; k < num_regions; ++k)
      h[k + 1] = h[k] + w_to_kw * power[k] / d_mdot;

    // Solve for pressure from outlet to inlet
    p[num_regions] = d_p_exit;
    for (int k = num_regions; k > 0; --k) {
      // Compute density at top and bottom of level from current
      // enthalpy and pressure
      double rho_high = iapws::rho_from_p_h(p[k], h[k]);
      double rho_low = iapws::rho_from_p_h(p[k - 1], h[k - 1]);

      // Compute velocities from densities
      double u_high = d_mdot / (rho_high * d_area);
      double u_low = d_mdot / (rho_low * d_area);

      // Compute new pressure
      constexpr double pa_to_mpa = 1e-6;
      p[k - 1] = p[k] + pa_to_mpa * ((d_mdot / d_area) * (u_high - u_low) +
                                     g * d_delta_z[k - 1] * rho_low);
    }

    // Check convergence on h and p
    h_diff = 0.0;
    p_diff = 0.0;
    double h_norm = 0.0;
    double p_norm = 0.0;
    for (int k = 0; k < num_regions + 1; ++k) {
      h_diff += std::pow(h[k] - h_old[k], 2);
      p_diff += std::pow(p[k] - p_old[k], 2);
      h_norm += std::pow(h[k], 2);
      p_norm += std::pow(p[k], 2);
    }
    h_diff = std::sqrt(h_diff) / std::sqrt(h_norm);
    p_diff = std::sqrt(p_diff) / std::sqrt(p_norm);

    converged = (h_diff < d_tol) && (p_diff < d_tol);
    if (converged)
      break;

    if (d_verbosity == HIGH) {
      std::cout << "At iteration " << it
                << ", relative enthalpy "
                   "difference is "
                << h_diff
                << " and relative pressure "
                   "difference is "
                << p_diff << std::endl;
    }
  }

  if (d_verbosity >= LOW) {
    if (converged) {
      std::cout << "Subchannel solver converged after " << it << " iterations"
                << std::endl;
    } else {
      std::cout << "Subchannel solver failed to converge after " << it << " iterations"
                << std::endl;
      std::cout << "Final relative enthalpy change is " << h_diff
                << " and relative pressure change is " << p_diff << std::endl;
    }
  }

  // Conversion from kg/m^3 to g/cm^3
  constexpr double kgm3_2_gcm3 = 0.001;

  // Compute temperature and density from enthalpy and pressure
  for (int k = 0; k < num_regions; ++k) {
    // Compute region-averaged quantities
    double h_mean = 0.5 * (h[k] + h[k + 1]);
    double p_mean = 0.5 * (p[k] + p[k + 1]);

    temperature[k] = iapws::T_from_p_h(p_mean, h_mean);
    density[k] = kgm3_2_gcm3 * iapws::rho_from_p_h(p_mean, h_mean);
  }
}

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
// end of Single_Pin_Subchannel.cpp
//---------------------------------------------------------------------------//
