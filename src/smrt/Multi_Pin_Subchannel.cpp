#include <algorithm>

#include "smrt/Multi_Pin_Subchannel.h"

#include "Nemesis/harness/Soft_Equivalence.hh"
#include <gsl/gsl>

namespace enrico {
//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
Multi_Pin_Subchannel::Multi_Pin_Subchannel(SP_Assembly assembly,
                                           RCP_PL parameters,
                                           const std::vector<double>& dz)
  : d_assembly(assembly)
{
  Expects(d_assembly != nullptr);

  double height = std::accumulate(dz.begin(), dz.end(), 0.0);
  Expects(nemesis::soft_equiv(height, d_assembly->height()));

  double mdot_per_area = parameters->get("mass_flow_rate", 0.4);

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

  auto inlet_temp = parameters->get("inlet_temperature", 565.0);
  auto exit_press = parameters->get("exit_pressure", 1.52e7);

  // Build single channel solver
  d_pin_subchannel = std::make_shared<Single_Pin_Subchannel>(parameters, dz);
  d_pin_subchannel->set_inlet_temperature(inlet_temp);
  d_pin_subchannel->set_exit_pressure(exit_press);

  // set up solution arrays
  generate_arrays();
}

void Multi_Pin_Subchannel::generate_arrays()
{
  int pins_x = d_assembly->num_pins_x();
  int pins_y = d_assembly->num_pins_y();

  pin_temps.resize(pins_x * pins_y * d_Nz);
  pin_densities.resize(pins_x * pins_y * d_Nz);
  pin_powers.resize(pins_x * pins_y * d_Nz);
}

//---------------------------------------------------------------------------//
// Solve subchannel equations over all pins
//---------------------------------------------------------------------------//
void Multi_Pin_Subchannel::solve(const std::vector<double>& powers)
{
  pin_powers = powers;

  // Convenience function to compute pin index
  auto pin_index = [pins_x, pins_y](int ix, int iy, int iz) {
    Expects(ix >= 0);
    Expects(iy >= 0);
    Expects(ix < pins_x);
    Expects(iy < pins_y);
    return ix + pins_x * (iy + pins_y * iz);
  };

  // Solve in each channel
  std::fill(pin_temps.begin(), pin_temps.end(), 0.0);
  std::fill(pin_densities.begin(), pin_densities.end(), 0.0);
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
          channel_power[iz] += pin_powers[pin_index(ix, iy, iz)];

        // Upper left
        if (ix > 0 && iy < pins_y)
          channel_power[iz] += pin_powers[pin_index(ix - 1, iy, iz)];

        // Lower right
        if (ix < pins_x && iy > 0)
          channel_power[iz] += pin_powers[pin_index(ix, iy - 1, iz)];

        // Lower left
        if (ix > 0 && iy > 0)
          channel_power[iz] += pin_powers[pin_index(ix - 1, iy - 1, iz)];

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
          pin_temps[idx] += 0.25 * channel_temp[iz];
          pin_densities[idx] += 0.25 * channel_density[iz];
        }

        // Upper left
        if (ix > 0 && iy < pins_y) {
          int idx = pin_index(ix - 1, iy, iz);
          pin_temps[idx] += 0.25 * channel_temp[iz];
          pin_densities[idx] += 0.25 * channel_density[iz];
        }

        // Lower right
        if (ix < pins_x && iy > 0) {
          int idx = pin_index(ix, iy - 1, iz);
          pin_temps[idx] += 0.25 * channel_temp[iz];
          pin_densities[idx] += 0.25 * channel_density[iz];
        }

        // Lower left
        if (ix > 0 && iy > 0) {
          int idx = pin_index(ix - 1, iy - 1, iz);
          pin_temps[idx] += 0.25 * channel_temp[iz];
          pin_densities[idx] += 0.25 * channel_density[iz];
        }
      }
    }
  }
}

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
// end of Multi_Pin_Subchannel.cpp
//---------------------------------------------------------------------------//
