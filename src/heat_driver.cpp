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

  xt::xtensor<double, 1> zs = xt::linspace(1, 5, 2);
  xt::xtensor<double, 1> rs = xt::linspace(5, 15, 2);
  VisualizationPin vpin(0.0, 0.0, zs, rs, 10);
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

  xt::xtensor<int, 4> cells = vpin.cells();
  int num_cells = cells.shape()[0]*cells.shape()[1]*cells.shape()[2];
  // REPLACE WITH ACCUMULATOR
  int num_entries = xt::where(xt::view(cells, xt::all(), xt::all(), xt::all(), xt::range(1, xt::placeholders::_)) >= 0).size();

  fh << "CELLS " << num_cells << " " << num_entries << "\n";

  int conn_size = cells.shape()[3] - 1;
  i = 0;
  for (auto c : xt::view(cells, xt::all(), xt::all(), xt::all(), xt::range(1, xt::placeholders::_))) {
    i++;

    // write cell value
    if ( c >= 0 ) { fh << c; }

    // start new row if starting a new element
    if (i % conn_size == 0) {
      fh << "\n";
    } else {
      if ( c >= 0) { fh << " "; }
    }
  }

  fh << "CELL_TYPES " << num_cells << "\n";
  for (auto v : xt::view(cells, xt::all(), xt::all(), xt::all(), 0)) {
    fh << v << "\n";
  }

  fh.close();

  return;
}

xt::xtensor<double, 3> SurrogateHeatDriver::VisualizationPin::points() {

  xt::xtensor<double, 3> pnts_out = xt::zeros<double>({axial_divs_ + 1, points_per_plane_, 3});

  xt::xtensor<double, 1> x = xt::zeros<double>({points_per_plane_});
  xt::xtensor<double, 1> y = xt::zeros<double>({points_per_plane_});

  xt::xtensor<double, 1> theta;
  theta = xt::linspace<double>(0., 2.*openmc::PI, t_res_ + 1);

  // first point is the pin center, start at one
  for(int i = 0; i < radial_divs_; i++) {
    double ring_rad = r_grid_(i);
    std::cout << "Ring radius: " << ring_rad << std::endl;
    for (int j = 1; j <= t_res_; j++) {
      int idx = i * t_res_ + j;
      x(idx) = ring_rad * std::cos(theta[j]);
      y(idx) = ring_rad * std::sin(theta[j]);
    }
  }

  // translate
  x += x_;
  y += y_;

  for (int i = 0; i < z_grid_.size(); i++) {
    double z = z_grid_[i];
    xt::view(pnts_out, i, xt::all(), 0) = x;
    xt::view(pnts_out, i, xt::all(), 1) = y;
    xt::view(pnts_out, i, xt::all(), 2) = z;
  }

  return pnts_out;
}

xt::xtensor<int, 4> SurrogateHeatDriver::VisualizationPin::cells() {
  int n_divs = z_grid_.size() - 1;
  int n_points = cells_per_plane_ + 1;
  xt::xtensor<int, 4> cells_out({axial_divs_, radial_divs_, t_res_, 10});

  // first ring is wedges
  xt::view(cells_out, xt::all(), 0, xt::all(), 0) = 13;
  xt::view(cells_out, xt::all(), 0, xt::all(), 1) = 6;

  // the rest are hexes
  xt::view(cells_out, xt::all(), xt::range(1, xt::placeholders::_), xt::all(), 0) = 12;
  xt::view(cells_out, xt::all(), xt::range(1, xt::placeholders::_), xt::all(), 1) = 8;

  xt::xtensor<int, 3> base = xt::zeros<int>({radial_divs_, t_res_, 8});

  // inner ring
  // cell connectivity for the first z level
  xt::view(base, 0, xt::all(), 1) = xt::arange(1, t_res_ + 1);
  xt::view(base, 0, xt::all(), 2) = xt::arange(2, t_res_ + 2);
  // adjust last cell
  xt::view(base, 0, t_res_ - 1, 2) = 1;
  xt::view(base, 0, xt::all(), 3) = n_points;
  xt::view(base, 0, xt::all(), 4) = xt::view(base, 0, xt::all(), 1) + points_per_plane_;
  xt::view(base, 0, xt::all(), 5) = xt::view(base, 0, xt::all(), 2) + points_per_plane_;


  xt::xtensor<int, 3> radial_base = xt::zeros<int>({1, t_res_, 8});
  xt::view(radial_base, xt::all(), xt::all(), 0) = xt::arange(1, t_res_ + 1);
  xt::view(radial_base, xt::all(), xt::all(), 1) = xt::arange(2, t_res_ + 2);
  xt::view(radial_base, xt::all(), xt::all(), 2) = xt::view(radial_base, xt::all(), xt::all(), 1) + t_res_;
  xt::view(radial_base, xt::all(), xt::all(), 3) = xt::view(radial_base, xt::all(), xt::all(), 0) + t_res_;
  xt::view(radial_base, xt::all(), t_res_ - 1, 1) = 1;
  xt::view(radial_base, xt::all(), t_res_ - 1, 2) = t_res_ + 1;

  xt::view(radial_base, xt::all(), xt::all(), 4) = xt::view(radial_base, xt::all(), xt::all(), 0) + points_per_plane_;
  xt::view(radial_base, xt::all(), xt::all(), 5) = xt::view(radial_base, xt::all(), xt::all(), 1) + points_per_plane_;
  xt::view(radial_base, xt::all(), xt::all(), 6) = xt::view(radial_base, xt::all(), xt::all(), 2) + points_per_plane_;
  xt::view(radial_base, xt::all(), xt::all(), 7) = xt::view(radial_base, xt::all(), xt::all(), 3) + points_per_plane_;

  // other rings
  xt::view(base, xt::range(1, xt::placeholders::_), xt::all(), xt::all()) = radial_base;
  for (int i = 1; i < radial_divs_; i++) {
    int start_idx = (i-1) * t_res_;
    xt::view(base, i, xt::all(), xt::all()) += start_idx;
  }


  for(int i = 0; i < z_grid_.size() - 1; i++) {
    xt::view(cells_out, i, xt::all(), xt::all(), xt::range(2, 10)) = base;
    // increment connectivity
    base += points_per_plane_;
  }

  // first ring should be wedges only, invalidate last two entries
  xt::view(cells_out, xt::all(), 0, xt::all(), xt::range(8,10)) = -1;

  return cells_out;
}

} // namespace stream
