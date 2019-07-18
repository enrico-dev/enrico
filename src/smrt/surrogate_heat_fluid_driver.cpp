#include <algorithm>

#include "enrico/utils.h"
#include "smrt/SurrogateHeatFluidDriver.h"

#include <gsl/gsl>

namespace enrico {
//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
SurrogateHeatFluidDriver::SurrogateHeatFluidDriver(SP_Assembly assembly,
                                                   RCP_PL subchannel_parameters,
                                                   RCP_PL conduction_parameters,
                                                   const std::vector<double>& dz)
  : d_assembly(assembly)
{
  Expects(d_assembly != nullptr);

  double height = std::accumulate(dz.begin(), dz.end(), 0.0);
  Expects(soft_equiv(height, d_assembly->height()));

  double mdot_per_area = subchannel_parameters->get("mass_flow_rate", 0.4);

  // Index extents
  d_Nx = d_assembly->num_pins_x() + 1;
  d_Ny = d_assembly->num_pins_y() + 1;
  d_Nz = dz.size();

  // Compute channel areas
  d_areas.resize(d_Nx * d_Ny, 0.0);

  // Loop over pins, not channels
  for (int iy = 0; iy < d_assembly->num_pins_y(); ++iy) {
    for (int ix = 0; ix < d_assembly->num_pins_x(); ++ix) {
      // Assign 1/4 of area to each neighboring channel
      double this_area = d_assembly->flow_area(ix, iy);

      d_areas[channel_index(ix, iy)] += 0.25 * this_area;
      d_areas[channel_index(ix + 1, iy)] += 0.25 * this_area;
      d_areas[channel_index(ix, iy + 1)] += 0.25 * this_area;
      d_areas[channel_index(ix + 1, iy + 1)] += 0.25 * this_area;
    }
  }

  // Compute mass flow rates in each channel
  d_mdots.resize(d_Nx * d_Ny, 0.0);
  for (int channel = 0; channel < d_Nx * d_Ny; ++channel)
    d_mdots[channel] = mdot_per_area * d_areas[channel];

  auto inlet_temp = subchannel_parameters->get("inlet_temperature", 565.0);
  pressure_bc_ = subchannel_parameters->get("exit_pressure", 15.2);

  // Build single channel solver
  d_pin_subchannel = std::make_unique<Single_Pin_Subchannel>(subchannel_parameters, dz);
  d_pin_subchannel->set_inlet_temperature(inlet_temp);
  d_pin_subchannel->set_exit_pressure(pressure_bc_);

  // Build single pin solver
  d_pin_conduction = std::make_unique<Single_Pin_Conduction>(conduction_parameters, dz);
  d_pin_conduction->set_fuel_radius(d_assembly->fuel_radius());
  d_pin_conduction->set_clad_radius(d_assembly->clad_radius());

  // set up solution arrays
  generate_arrays();
}

void SurrogateHeatFluidDriver::generate_arrays()
{
  int pins_x = d_assembly->num_pins_x();
  int pins_y = d_assembly->num_pins_y();

  d_pin_temps.resize(pins_x * pins_y * d_Nz);
  d_pin_densities.resize(pins_x * pins_y * d_Nz);
  d_pin_powers.resize(pins_x * pins_y * d_Nz);
}

void SurrogateHeatFluidDriver::solve(const std::vector<double>& powers)
{
  solve_fluid(powers);
  solve_heat(powers, d_pin_temps, d_solid_temps);
}

//---------------------------------------------------------------------------//
// Solve subchannel equations over all pins
//---------------------------------------------------------------------------//
void SurrogateHeatFluidDriver::solve_fluid(const std::vector<double>& powers)
{
  d_pin_powers = powers;

  // Convenience function to compute pin index
  auto pin_index = [pins_x, pins_y](int ix, int iy, int iz) {
    Expects(ix >= 0);
    Expects(iy >= 0);
    Expects(ix < pins_x);
    Expects(iy < pins_y);
    return ix + pins_x * (iy + pins_y * iz);
  };

  // Solve in each channel
  std::fill(d_pin_temps.begin(), d_pin_temps.end(), 0.0);
  std::fill(d_pin_densities.begin(), d_pin_densities.end(), 0.0);
  std::vector<double> channel_power(d_Nz);
  std::vector<double> channel_temp(d_Nz);
  std::vector<double> channel_density(d_Nz);
  for (int iy = 0; iy < d_Ny; ++iy) {
    for (int ix = 0; ix < d_Nx; ++ix) {
      // Get power for this channel
      for (int iz = 0; iz < d_Nz; ++iz) {
        channel_power[iz] = 0.0;

        // Upper right
        if (ix < pins_x && iy < pins_y)
          channel_power[iz] += d_pin_powers[pin_index(ix, iy, iz)];

        // Upper left
        if (ix > 0 && iy < pins_y)
          channel_power[iz] += d_pin_powers[pin_index(ix - 1, iy, iz)];

        // Lower right
        if (ix < pins_x && iy > 0)
          channel_power[iz] += d_pin_powers[pin_index(ix, iy - 1, iz)];

        // Lower left
        if (ix > 0 && iy > 0)
          channel_power[iz] += d_pin_powers[pin_index(ix - 1, iy - 1, iz)];

        channel_power[iz] *= 0.25;
      }

      // Solve for this channel
      int channel_id = channel_index(ix, iy);
      d_pin_subchannel->set_channel_area(d_areas[channel_id]);
      d_pin_subchannel->set_mass_flow_rate(d_mdots[channel_id]);
      d_pin_subchannel->solve(channel_power, channel_temp, channel_density);

      // Add contribution of temp and density into pin-centered values
      for (int iz = 0; iz < d_Nz; ++iz) {
        // Upper right
        if (ix < pins_x && iy < pins_y) {
          int idx = pin_index(ix, iy, iz);
          d_pin_temps[idx] += 0.25 * channel_temp[iz];
          d_pin_densities[idx] += 0.25 * channel_density[iz];
        }

        // Upper left
        if (ix > 0 && iy < pins_y) {
          int idx = pin_index(ix - 1, iy, iz);
          d_pin_temps[idx] += 0.25 * channel_temp[iz];
          d_pin_densities[idx] += 0.25 * channel_density[iz];
        }

        // Lower right
        if (ix < pins_x && iy > 0) {
          int idx = pin_index(ix, iy - 1, iz);
          d_pin_temps[idx] += 0.25 * channel_temp[iz];
          d_pin_densities[idx] += 0.25 * channel_density[iz];
        }

        // Lower left
        if (ix > 0 && iy > 0) {
          int idx = pin_index(ix - 1, iy - 1, iz);
          d_pin_temps[idx] += 0.25 * channel_temp[iz];
          d_pin_densities[idx] += 0.25 * channel_density[iz];
        }
      }
    }
  }
}

void SurrogateHeatFluidDriver::solve_heat(const std::vector<double>& power,
                                          const std::vector<double>& channel_temp,
                                          std::vector<double>& fuel_temp)
{
  int Nx = d_assembly->num_pins_x();
  int Ny = d_assembly->num_pins_y();
  int N = Nx * Ny * d_Nz;

  Expects(power.size() == N);
  Expects(channel_temp.size() == N);
  Expects(fuel_temp.size() == N);

  // Storage for single pin values
  std::vector<double> pin_power(d_Nz);
  std::vector<double> pin_channel_temp(d_Nz);
  std::vector<double> pin_fuel_temp(d_Nz);

  // Loop over pins
  for (int iy = 0; iy < Ny; ++iy) {
    for (int ix = 0; ix < Nx; ++ix) {
      // Copy assembly data into single-pin containers
      for (int iz = 0; iz < d_Nz; ++iz) {
        int assembly_idx = ix + Nx * (iy + Ny * iz);
        pin_power[iz] = power[assembly_idx];
        pin_channel_temp[iz] = channel_temp[assembly_idx];
      }

      if (d_assembly->pin_type(ix, iy) == Assembly_Model::FUEL) {
        // Solve conduction in pin
        d_pin_conduction->solve(pin_power, pin_channel_temp, pin_fuel_temp);
      } else {
        for (int iz = 0; iz < d_Nz; ++iz)
          Expects(pin_power[iz] == 0);

        // "Fuel" temp is same as channel temp in guide tubes
        std::copy(
          pin_channel_temp.begin(), pin_channel_temp.end(), pin_fuel_temp.begin());
      }

      // Copy fuel temp data back to assembly container
      for (int iz = 0; iz < d_Nz; ++iz) {
        int assembly_idx = ix + Nx * (iy + Ny * iz);
        fuel_temp[assembly_idx] = pin_fuel_temp[iz];
      }
    }
  }
}

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
// end of SurrogateHeatFluidDriver.cpp
//---------------------------------------------------------------------------//
