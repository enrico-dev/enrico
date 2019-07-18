#include <cmath>
#include <iostream>

#include "enrico/utils.h"
#include "smrt/Single_Pin_Conduction.h"

#include "Nemesis/utils/String_Functions.hh"
#include "openmc/constants.h"

namespace enrico {
//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
Single_Pin_Conduction::Single_Pin_Conduction(RCP_PL& parameters,
                                             const std::vector<double>& delta_z)
  : d_fuel_radius(0.0)
  , d_clad_radius(0.0)
  , d_delta_z(delta_z)
{
  // Read parameters from parameter list
  d_delta_r_fuel = parameters->get("delta_r_fuel", 0.05);
  d_delta_r_clad = parameters->get("delta_r_clad", 0.02);
  d_k_fuel = parameters->get("fuel_conductivity", 0.0287);
  d_k_clad = parameters->get("clad_conductivity", 0.215);
}

//---------------------------------------------------------------------------//
// Solve for fuel temperature
//---------------------------------------------------------------------------//
void Single_Pin_Conduction::solve(const std::vector<double>& power,
                                  const std::vector<double>& channel_temp,
                                  std::vector<double>& fuel_temp)
{
  Expects(power.size() == d_delta_z.size());
  Expects(channel_temp.size() == d_delta_z.size());
  Expects(fuel_temp.size() == d_delta_z.size());
  Expects(d_fuel_radius > 0.0);
  Expects(d_clad_radius > d_fuel_radius);

  using openmc::PI;

  int num_levels = d_delta_z.size();

  // Number of radial rings in each region
  int num_rings_fuel = std::ceil(d_fuel_radius / d_delta_r_fuel);
  int num_rings_clad = std::ceil((d_clad_radius - d_fuel_radius) / d_delta_r_clad);
  int num_rings = num_rings_fuel + num_rings_clad;
  Expects(num_rings > 2);

  // Compute delta r in fuel and clad
  std::vector<double> radius(num_rings + 1);
  double dr = d_fuel_radius / static_cast<double>(num_rings_fuel);
  radius[0] = 0.0;
  for (int r = 0; r < num_rings_fuel; ++r)
    radius[r + 1] = radius[r] + dr;
  dr = (d_clad_radius - d_fuel_radius) / static_cast<double>(num_rings_clad);
  for (int r = num_rings_fuel; r < num_rings; ++r)
    radius[r + 1] = radius[r] + dr;
  Expects(soft_equiv(radius[num_rings], d_clad_radius));

  // Assign thermal conductivity in fuel and clad
  std::vector<double> k(num_rings);
  std::fill(k.begin(), k.begin() + num_rings_fuel, d_k_fuel);
  std::fill(k.begin() + num_rings_fuel, k.end(), d_k_clad);

  // Compute area
  std::vector<double> area(num_rings);
  for (int r = 0; r < num_rings; ++r) {
    area[r] = PI * (radius[r + 1] * radius[r + 1] - radius[r] * radius[r]);
  }
  // Check areas
  double fuel_area = std::accumulate(area.begin(), area.begin() + num_rings_fuel, 0.0);
  Expects(soft_equiv(fuel_area, PI * d_fuel_radius * d_fuel_radius));

  double total_area = std::accumulate(area.begin(), area.end(), 0.0);
  Expects(soft_equiv(total_area, PI * d_clad_radius * d_clad_radius));

  auto rcp_matrix = Teuchos::rcp(new Matrix(num_rings, num_rings));
  auto rcp_rhs = Teuchos::rcp(new Vector(num_rings));
  auto rcp_soln = Teuchos::rcp(new Vector(num_rings));
  Solver solver;
  double dr_lo, dr_hi;
  for (int l = 0; l < num_levels; ++l) {
    // Dereference RCPs to get objects
    auto& matrix = *rcp_matrix;
    auto& rhs = *rcp_rhs;
    auto& soln = *rcp_soln;

    matrix.putScalar(0.0);
    rhs.putScalar(0.0);
    soln.putScalar(0.0);

    // First cell
    dr_hi = 2 * radius[2] / radius[1];
    matrix(0, 0) = 2 * PI * k[0] * dr_hi;
    matrix(0, 1) = -2 * PI * k[0] * dr_hi;

    // Inner cells
    for (int r = 1; r < num_rings - 1; ++r) {
      dr_hi = 0.5 * (radius[r + 2] - radius[r]);
      dr_lo = 0.5 * (radius[r + 1] - radius[r - 1]);
      matrix(r, r - 1) = -2 * PI * k[r] * radius[r] / dr_lo;
      matrix(r, r) = 2 * PI * k[r] * (radius[r] / dr_lo + radius[r + 1] / dr_hi);
      matrix(r, r + 1) = -2 * PI * k[r] * radius[r + 1] / dr_hi;
    }

    // Last cell
    int r = num_rings - 1;
    dr_lo = 0.5 * (radius[r + 1] - radius[r - 1]);
    dr_hi = 0.5 * (radius[r + 1] - radius[r]);
    matrix(r, r - 1) = -2 * PI * k[r] * radius[r] / dr_lo;
    matrix(r, r) = 2 * PI * k[r] * (radius[r] / dr_lo + radius[r + 1] / dr_hi);

    // Right hand side
    for (int r = 0; r < num_rings_fuel; ++r)
      rhs[r] = power[l] * area[r] / d_delta_z[l] / fuel_area;

    rhs[num_rings - 1] += 2 * PI * k.back() * channel_temp[l] * d_clad_radius / dr_hi;

    // Solve linear system
    solver.setMatrix(rcp_matrix);
    solver.setVectors(rcp_soln, rcp_rhs);
    solver.equilibrateMatrix();
    solver.equilibrateRHS();
    int err = solver.solve();
    solver.unequilibrateLHS();
    Expects(err == 0);

    // Compute average fuel temperature for this level
    double ft = 0.0;
    for (int r = 0; r < num_rings_fuel; ++r) {
      ft += soln[r] * area[r];
    }
    fuel_temp[l] = ft / fuel_area;
  }
}

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
// end of Single_Pin_Conduction.cpp
//---------------------------------------------------------------------------//
