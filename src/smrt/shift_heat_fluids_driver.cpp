//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   shift_heat_fluids_driver.cpp
 * \author Steven Hamilton
 * \date   Fri Aug 10 09:00:11 2018
 * \brief  ShiftHeatFluidsDriver class definitions.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "smrt/shift_heat_fluids_driver.hpp"

// SCALE includes
#include "Nemesis/comm/Logger.hh"
#include "Nemesis/comm/global.hh"
#include "Nemesis/harness/Soft_Equivalence.hh"
#include "Omnibus/config.h"

// vendored includes
#include <gsl/gsl>

// enrico includes
#include "smrt/shift_driver.h"

namespace enrico {
//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
ShiftHeatFluidsDriver::ShiftHeatFluidsDriver(SP_Assembly assembly,
                                             RCP_PL params,
                                             const Vec_Dbl& z_edges)
  : d_assembly(assembly)
{
  Expects(assembly != nullptr);

  Expects(nemesis::soft_equiv(z_edges.back(), d_assembly->height()));

  Vec_Dbl dz(z_edges.size() - 1);
  for (int edge = 0; edge < dz.size(); ++edge)
    dz[edge] = z_edges[edge + 1] - z_edges[edge];

  // Resize storage for solutions
  d_Nx = d_assembly->num_pins_x();
  d_Ny = d_assembly->num_pins_y();
  d_Nz = dz.size();
  d_N = d_Nx * d_Ny * d_Nz;

  d_power.resize(d_N);
  d_fuel_temperature.resize(d_N);
  d_coolant_temperature.resize(d_N);
  d_coolant_density.resize(d_N);

  d_power_magnitude = params->get("total_power", 70000.0 * d_Nx * d_Ny);
  d_tol = params->get("tolerance", 1e-4);
  d_max_iters = params->get("max_iters", 25);
  d_damping = params->get("relaxation_factor", 1.0);
  d_verbosity = params->get("verbosity", std::string("low"));

  // Build fluids-heat solver
  auto subchannel_params = Teuchos::sublist(params, "Subchannel");
  auto conduction_params = Teuchos::sublist(params, "Conduction");
  d_heat_fluid = std::make_shared<SurrogateHeatFluidDriver>(
    assembly, subchannel_params, conduction_params, dz);

  // Build neutronics solver (surrogate diffusion or Shift MC)
  auto neutronics_params = Teuchos::sublist(params, "Neutronics");
  auto shift_input = neutronics_params->get<std::string>("shift_input");
  d_neutronics = std::make_shared<Shift_Driver>(assembly, shift_input, z_edges);
}

//---------------------------------------------------------------------------//
// Solve
//---------------------------------------------------------------------------//
void ShiftHeatFluidsDriver::solve()
{
  // Set initial guess for power
  std::fill(d_power.begin(), d_power.end(), 0.0);
  for (int jpin = 0; jpin < d_assembly->num_pins_y(); ++jpin) {
    for (int ipin = 0; ipin < d_assembly->num_pins_x(); ++ipin) {
      if (d_assembly->pin_type(ipin, jpin) == Assembly::FUEL) {
        for (int k = 0; k < d_Nz; ++k) {
          d_power[d_assembly->index(ipin, jpin, k)] = 1.0;
        }
      }
    }
  }
  double power_nrm = std::accumulate(d_power.begin(), d_power.end(), 0.0);
  for (auto& val : d_power)
    val *= d_power_magnitude / power_nrm;

  // Previous iteration power for convergence checking
  Vec_Dbl old_power;

  for (int it = 0; it < d_max_iters; ++it) {
    old_power = d_power;

    // solve the fluid and solid T/H equations
    d_heat_fluid->solve(d_power);

    // get the needed portions of the T/H solution
    d_coolant_temperature = d_heat_fluid->fluid_temperature();
    d_coolant_density = d_heat_fluid->density();
    d_fuel_temperature = d_heat_fluid->solid_temperature();

    d_neutronics->solve(d_fuel_temperature, d_coolant_density, d_power);

    // Scale power
    double power_nrm = 0.0;
    for (const auto& val : d_power)
      power_nrm += val;
    for (auto& val : d_power)
      val *= d_power_magnitude / power_nrm;

    // Compute relative difference in power
    double power_diff = 0.0;
    for (int i = 0; i < d_N; ++i)
      power_diff += std::abs(d_power[i] - old_power[i]);
    power_diff /= d_power_magnitude;

    if (nemesis::node() == 0) {
      std::cout << "Relative power difference at iteration " << it << " = " << power_diff
                << std::endl;
    }

    // Print solutions if requested
    if (d_verbosity == "high" && nemesis::node() == 0) {
      std::cout << std::endl;
      std::cout << "Coolant density: ";
      for (auto val : d_coolant_density)
        std::cout << val << " ";

      std::cout << std::endl << std::endl;
      std::cout << "Coolant temp: ";
      for (auto val : d_coolant_temperature)
        std::cout << val << " ";

      std::cout << std::endl << std::endl;
      std::cout << "Fuel temp: ";
      for (auto val : d_fuel_temperature)
        std::cout << val << " ";

      std::cout << std::endl << std::endl;
      std::cout << "Power : ";
      for (auto val : d_power)
        std::cout << val << " ";
      std::cout << std::endl << std::endl;
    }

    // Check for convergence
    if (power_diff / d_damping < d_tol) {
      if (nemesis::node() == 0) {
        std::cout << "Coupled physics solve converged after " << it << " iterations"
                  << std::endl;
      }
      break;
    }

    // Apply relaxation
    for (int i = 0; i < d_N; ++i) {
      d_power[i] = d_damping * d_power[i] + (1.0 - d_damping) * old_power[i];
    }
  }
}

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
// end of smrt/shift_heat_fluids_driver.cpp
//---------------------------------------------------------------------------//
