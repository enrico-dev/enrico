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

using namespace xt::placeholders;

xt::xtensor<double, 2> create_ring_points(double radius,
                                          int t_resolution) {
  xt::xtensor<double, 1> x = xt::zeros<double>({t_resolution});
  xt::xtensor<double, 1> y = xt::zeros<double>({t_resolution});

  xt::xtensor<double, 1> theta;
  theta = xt::linspace<double>(0., 2.*openmc::PI, t_resolution + 1);

  x = radius * xt::cos(theta);
  y = radius * xt::sin(theta);

  xt::xtensor<double, 2> out = xt::zeros<double>({2, t_resolution});
  xt::view(out, 0, xt::all()) = x;
  xt::view(out, 1, xt::all()) = x;

  return out;
}

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

void SurrogateHeatDriver::to_vtk(std::string filename)
{
  std::cout << "Writing VTK file: " << filename << "...\n";

  xt::xtensor<double, 1> zs = xt::linspace(1, 5, 5);
  xt::xtensor<double, 1> rs = xt::linspace(5, 15, 5);

  // create a pin
  int radial_resolution = 50;
  VisualizationPin vpin(pin_centers_(0,0),
                        pin_centers_(0,1),
                        z_,
                        r_grid_fuel_,
                        radial_resolution);
  xt::xtensor<double, 3> pin_points = vpin.points();

  // open vtk file
  std::ofstream fh(filename, std::ofstream::out);

  // write header
  fh << "# vtk DataFile Version 2.0\n";
  fh << "No comment\nASCII\nDATASET UNSTRUCTURED_GRID\n";
  fh << "POINTS " << vpin.num_points() << " float\n";
  xt::xtensor<double, 1> points_flat = xt::flatten(pin_points);
  for (auto p = points_flat.begin(); p != points_flat.end(); p+=3) {
    fh << *p << " " << *(p+1) << " " << *(p+2) << "\n";
  }

  // generate cell connectivity
  xt::xtensor<int, 4> cells = vpin.cells();

  // separate cell types and cell entries
  xt::xtensor<int, 4> cell_types = xt::view(cells, xt::all(), xt::all(), xt::all(), xt::range(0,1));
  cells = xt::view(cells, xt::all(), xt::all(), xt::all(), xt::range(1, _));

  fh << "\nCELLS " << vpin.num_cells() << " " << vpin.num_entries() << "\n";

  // write points to file
  xt::xtensor<int, 1> cells_flat = xt::flatten(cells);
  int conn_size = vpin.conn_entry_size();
  for (auto c = cells_flat.begin(); c != cells_flat.end(); c += conn_size) {
    for (int i = 0; i < conn_size; i++) {
      auto val = *(c+i);
      if (val >= 0) { fh << val << " "; }
    }
    fh << "\n";
  }

  // write cell types
  fh << "\nCELL_TYPES " << vpin.num_cells() << "\n";
  for (auto v : cell_types) {
    fh << v << "\n";
  }

  // data header
  fh << "CELL_DATA " << vpin.num_cells() << "\n";

  // temperature data
  fh << "SCALARS TEMPERATURE double 1\n";
  fh << "LOOKUP_TABLE default\n";
  for (int i = 0; i < temperature_.shape()[1]; i++) {
    for (int j = 0; j < vpin.radial_divs_; j++) {
      for (int k = 0; k < radial_resolution; k++) {
        fh << temperature_(0, i, j) << "\n";
      }
    }
  }

  // fission source data
  fh << "SCALARS SOURCE double 1\n";
  fh << "LOOKUP_TABLE default\n";
  for (int i = 0; i < source_.shape()[1]; i++) {
    for (int j = 0; j < vpin.radial_divs_; j++) {
      for (int k = 0; k < radial_resolution; k++) {
        fh << source_(0, i, j) << "\n";
      }
    }
  }

  fh.close();

  return;
}

SurrogateHeatDriver::VisualizationPin::VisualizationPin(double x, double y,
                                   xt::xtensor<double, 1> z_grid,
                                   xt::xtensor<double, 1> r_grid,
                                   int t_res)
    : x_(x), y_(y), z_grid_(z_grid), r_grid_(r_grid), t_res_(t_res)
    {
      // adjust to remove extra radial ring if needed
      if (r_grid[0] == 0) {
        r_grid_ = xt::view(r_grid_, xt::range(1,_));
      }
      axial_divs_ = z_grid_.size() - 1;
      radial_divs_ = r_grid_.size();
      points_per_plane_ = radial_divs_ * t_res_ + 1;
      cells_per_plane_ = points_per_plane_ - 1;
    }

xt::xtensor<double, 3> SurrogateHeatDriver::VisualizationPin::points() {

  xt::xarray<double> pnts_out = xt::zeros<double>({axial_divs_ + 1,
                                                  points_per_plane_,
                                                  3});

  xt::xtensor<double, 1> x = xt::zeros<double>({points_per_plane_});
  xt::xtensor<double, 1> y = xt::zeros<double>({points_per_plane_});

  xt::xtensor<double, 1> theta;
  theta = xt::linspace<double>(0., 2.*openmc::PI, t_res_ + 1);

  // first point is the pin center, start at one
  for(int i = 0; i < radial_divs_; i++) {
    double ring_rad = r_grid_(i);
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
    xt::view(pnts_out, i, xt::all(), 0) = x;
    xt::view(pnts_out, i, xt::all(), 1) = y;
    xt::view(pnts_out, i, xt::all(), 2) = z_grid_[i];
  }

  return pnts_out;
}

xt::xtensor<int, 4> SurrogateHeatDriver::VisualizationPin::cells() {
  // size output array
  xt::xtensor<int, 4> cells_out = xt::zeros<int>({axial_divs_,
                                                  radial_divs_,
                                                  t_res_,
                                                  HEX_SIZE_ + 2});

  // generate a base layer to be extended in Z
  xt::xtensor<int, 3> base = xt::zeros<int>({radial_divs_, t_res_, HEX_SIZE_});

  /// INNER RING \\\

  xt::xtensor<int, 2> inner_base = xt::zeros<int>({t_res_, 8});

  // cell connectivity for the first z level
  xt::view(inner_base, xt::all(), 1) = xt::arange(1, t_res_ + 1);
  xt::view(inner_base, xt::all(), 2) = xt::arange(2, t_res_ + 2);
  // adjust last cell for perioic condition
   xt::view(inner_base, t_res_ - 1, 2) = 1;

  // copy connectivity of first layer to the second
  xt::view(inner_base, xt::all(), xt::range(3,6)) =
    xt::view(inner_base, xt::all(), xt::range(0,3));
  // shift connectivity down one layer
  xt::strided_view(inner_base, {xt::all(), xt::range(3,6)}) += points_per_plane_;

  // set the inner_base
   xt::view(base, 0, xt::all(), xt::all()) = inner_base;

  /// OUTER RINGS \\\

  xt::xtensor<int, 2> radial_base = xt::zeros<int>({t_res_, 8});

  // setup connectivity of the first layer
  xt::view(radial_base, xt::all(), 0) = xt::arange(1, t_res_ + 1);
  xt::view(radial_base, xt::all(), 1) = xt::arange(2, t_res_ + 2);
  xt::view(radial_base, xt::all(), 2) = xt::view(radial_base, xt::all(), 1) + t_res_;
  xt::view(radial_base, xt::all(), 3) = xt::view(radial_base, xt::all(), 0) + t_res_;
  xt::view(radial_base, t_res_ - 1, 1) = 1;
  xt::view(radial_base, t_res_ - 1, 2) = t_res_ + 1;

  // copy connectivity of the first layer
  xt::view(radial_base, xt::all(), xt::range(4,8)) +=
    xt::view(radial_base, xt::all(), xt::range(0,4));
  // shift connectivity down one layer
  xt::view(radial_base, xt::all(), xt::range(4,8)) += points_per_plane_;

  // other rings
  xt::strided_view(base, {xt::range(1, _), xt::ellipsis()}) = radial_base;
  for (int i = 1; i < radial_divs_; i++) {
    // extend based on the starting index of the ring
    int start_idx = (i-1) * t_res_;
    xt::view(base, i, xt::all(), xt::all()) += start_idx;
  }

  // the inner ring is wedges
  xt::view(cells_out, xt::all(), 0, xt::all(), 0) = WEDGE_TYPE_;
  xt::view(cells_out, xt::all(), 0, xt::all(), 1) = WEDGE_SIZE_;

  // the rest are hexes
  xt::view(cells_out, xt::all(), xt::range(1, _), xt::all(), 0) = HEX_TYPE_;
  xt::view(cells_out, xt::all(), xt::range(1, _), xt::all(), 1) = HEX_SIZE_;

  // set all axial divs using base
  for(int i = 0; i < axial_divs_; i++) {
    // set layer and increment connectivity by number of points in axial div
    xt::view(cells_out, i, xt::all(), xt::all(), xt::range(2, 10)) = base;
    base += points_per_plane_;
  }

  // first ring should be wedges only, invalidate last two entries
  xt::view(cells_out, xt::all(), 0, xt::all(), xt::range(8,10)) = INVALID_CONN_;

  return cells_out;
}

} // namespace stream
