#include "stream/heat_driver.h"

#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xview.hpp"
#include "heat_xfer_backend.h"
#include "openmc/xml_interface.h"
#include "openmc/constants.h"

#include <iostream>
#include <fstream>

namespace stream {

SurrogateHeatDriver::SurrogateHeatDriver(MPI_Comm comm, pugi::xml_node node)
  : HeatFluidsDriver{comm}
{
  // Determine heat transfer solver parameters
  clad_inner_radius_ = node.child("clad_inner_radius").text().as_double();
  clad_outer_radius_ = node.child("clad_outer_radius").text().as_double();
  pellet_radius_ = node.child("pellet_radius").text().as_double();
  n_fuel_rings_ = node.child("fuel_rings").text().as_int();
  n_clad_rings_ = node.child("clad_rings").text().as_int();

  // Get pin locations
  // TODO: Switch to get_node_xarray on OpenMC update
  auto pin_locations = openmc::get_node_array<double>(node, "pin_centers");
  if (pin_locations.size() % 2 != 0) {
    throw std::runtime_error{"Length of <pin_centers> must be a multiple of two"};
  }

  // Convert to xtensor
  n_pins_ = pin_locations.size() / 2;
  std::vector<std::size_t> shape {n_pins_, 2};
  pin_centers_ = xt::adapt(pin_locations, shape);

  // Get z values
  // TODO: Switch to get_node_xarray on OpenMC update
  auto z_values = openmc::get_node_array<double>(node, "z");
  z_ = xt::adapt(z_values);
  n_axial_ = z_.size() - 1;

  // Heat equation solver tolerance
  tol_ = node.child("tolerance").text().as_double();

  // Initialize heat transfer solver
  generate_arrays();
};

void SurrogateHeatDriver::generate_arrays()
{
  // Make a radial grid for the fuel with equal spacing.
  r_grid_clad_ = xt::linspace<double>(
    clad_inner_radius_,
    clad_outer_radius_,
    n_clad_rings_ + 1
  );

  // Make a radial grid for the clad with equal spacing.
  r_grid_fuel_ = xt::linspace<double>(
    0, pellet_radius_, n_fuel_rings_ + 1);

  // Create empty arrays for source term and temperature
  source_ = xt::empty<double>({n_pins_, n_axial_, n_rings()});
  temperature_ = xt::empty<double>({n_pins_, n_axial_, n_rings()});
}

void SurrogateHeatDriver::solve_step()
{
  std::cout << "Solving heat equation...\n";
  // NuScale inlet temperature
  double T_co = 523.15;

  // Set initial temperature
  for (auto& T : temperature_) {
    T = T_co;
  }

  // Convert source to [W/m^3] as expected by Magnolia
  xt::xtensor<double, 3> q = 1e6*source_;

  // Convert radial grids to [m] as expected by Magnolia
  xt::xtensor<double, 1> r_fuel = 0.01*r_grid_fuel_;
  xt::xtensor<double, 1> r_clad = 0.01*r_grid_clad_;

  for (int i = 0; i < n_pins_; ++i) {
    for (int j = 0; j < n_axial_; ++j) {
      solve_steady_nonlin(&q(i, j, 0), T_co, r_fuel.data(), r_clad.data(),
        n_fuel_rings_, n_clad_rings_, tol_, &temperature_(i, j, 0));
    }
  }
}

void SurrogateHeatDriver::to_vtk(std::string filename,
                                 VTKData output_data)
{

  std::cout << "Writing VTK file: " << filename << "...\n";

  std::vector<double> zs = {0.0, 1.0, 2.0};
  VisualizationPin vpin(10.0, 10.0, 5.0, zs, 10);
  xt::xtensor<double, 3> pin_points = vpin.points();

  // open vtk file
  std::ofstream fh(filename, std::ofstream::out);

  // write header
  fh << "# vtk DataFile Version 2.0\n";
  fh << "No comment\nASCII\nDATASET UNSTRUCTURED_GRID\n";
  fh << "POINTS " << pin_points.shape()[0] * pin_points.shape()[1] << " float\n";

  int i = 0;
  for (auto p : pin_points) {
    i++;
    fh << p;
    if (i % 3 == 0) {
      fh << "\n";
    } else {
      fh <<  " ";
    }
  }

  xt::xtensor<int, 3> cells = vpin.cells();
  int num_cells = cells.shape()[0]*cells.shape()[1];
  fh << "CELLS " << num_cells << " " << num_cells * 7 << "\n";

  i = 0;
  for (auto c : cells) {
    i++;
    fh << c;
    if (i % 7 == 0) {
      fh << "\n";
    } else {
      fh << " ";
    }
  }

  fh << "CELL_TYPES " << num_cells << "\n";
  for (int j = 0; j < num_cells; j++) {
    fh << 13 << "\n";
  }

  fh.close();

  return;
}

xt::xtensor<double, 3> SurrogateHeatDriver::VisualizationPin::points() {
  int n_divs = z_grid.size() - 1;
  int points_per_plane = t_resolution + 1;

  xt::xtensor<double, 3> pnts_out = xt::zeros<double>({n_divs + 1, points_per_plane, 3});

  xt::xtensor<double, 1> x = xt::zeros<double>({points_per_plane});
  xt::xtensor<double, 1> y = xt::zeros<double>({points_per_plane});

  xt::xtensor<double, 1> theta;
  theta = xt::linspace<double>(0., 2.*openmc::PI, t_resolution + 1);

  // populate x and y values
  x(0, 0, 0) = 0;
  y(0, 0, 0) = 0;

  for (int i = 1; i < points_per_plane; i++) {
    x[i] = std::cos(theta[i]);
    y[i] = std::sin(theta[i]);
  }
  // scale
  x *= pin_radius;
  y *= pin_radius;
  // translate
  x += x_;
  y += y_;

  for (int i = 0; i < z_grid.size() ; i ++) {
    double z = z_grid[i];
    xt::view(pnts_out, i, xt::all(), 0) = x;
    xt::view(pnts_out, i, xt::all(), 1) = y;
    xt::view(pnts_out, i, xt::all(), 2) = z;
  }

  return pnts_out;
}

xt::xtensor<int, 3> SurrogateHeatDriver::VisualizationPin::cells() {
  int n_divs = z_grid.size() - 1;
  int n_points = t_resolution + 1;
  xt::xtensor<int, 3> cells_out({n_divs, t_resolution, 7});

  xt::view(cells_out, xt::all(), xt::all(), 0) = 6;

  xt::xtensor<int, 2> base = xt::zeros<int>({t_resolution, 6});

  // cell connectivity for the first z level
  xt::view(base, xt::all(), 1) = xt::arange(1, n_points);
  xt::view(base, xt::all(), 2) = xt::arange(2, n_points + 1);
  // adjust last cell
  xt::view(base, t_resolution-1, 2) = 1;
  xt::view(base, xt::all(), 3) = n_points;
  xt::view(base, xt::all(), 4) = xt::view(base, xt::all(), 1) + n_points;
  xt::view(base, xt::all(), 5) = xt::view(base, xt::all(), 2) + n_points;

  for(int i = 0; i < z_grid.size() - 1; i++) {
    base += n_points * i;
    xt::view(cells_out, i, xt::all(), xt::range(1, 7)) = base;
  }

  return cells_out;
}

} // namespace stream
