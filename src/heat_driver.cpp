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
  // create a pin
  int radial_resolution = 40;

  int n_axial_sections = z_.size() - 1;
  int n_radial_fuel_sections = r_grid_fuel_.size() - 1;
  int n_radial_clad_sections = r_grid_clad_.size() - 1;
  int n_radial_sections = n_radial_fuel_sections + n_radial_clad_sections;
  int n_sections_per_plane = n_radial_sections * radial_resolution;
  // Calculate some necessary values ahead of time
  // fuel points
  int fuel_points_per_plane  = n_radial_fuel_sections * radial_resolution + 1;
  int fuel_points = fuel_points_per_plane * z_.size();
  // cladding points
  int clad_points_per_plane = r_grid_clad_.size() * radial_resolution;
  int points_per_plane = fuel_points_per_plane + clad_points_per_plane;
  // total number of points in pin
  int total_points = points_per_plane * z_.size();

  // fuel elements, entries
  int num_mesh_elements = n_sections_per_plane * n_axial_sections;
  // wedge regions
  int num_fuel_entries_per_plane = radial_resolution * (WEDGE_SIZE_ + 1);
  // other radial regions
  num_fuel_entries_per_plane += radial_resolution * (HEX_SIZE_ + 1) * (n_radial_fuel_sections - 1);
  int num_clad_entries_per_plane = radial_resolution * (HEX_SIZE_ + 1) * n_radial_clad_sections;
  int num_entries_per_plane = num_fuel_entries_per_plane + num_clad_entries_per_plane;
  // cladding elements
  int num_entries = num_entries_per_plane * n_axial_sections;

  VisualizationPin vpin(pin_centers_(0,0),
                        pin_centers_(0,1),
                        z_,
                        r_grid_fuel_,
                        r_grid_clad_,
                        radial_resolution);


  // generate fuel and cladding points
  xt::xtensor<double, 3> pin_points = vpin.fuel_points();
  xt::xtensor<double, 3> clad_points = vpin.clad_points();

  // generate mesh element connectivity for fuel
  xt::xtensor<int, 4> cells = vpin.fuel_connectivity();
  xt::xtensor<int, 4> cell_types = xt::view(cells, xt::all(), xt::all(), xt::all(), xt::range(0,1));
  cells = xt::view(cells, xt::all(), xt::all(), xt::all(), xt::range(1, _));

  // generate mesh element connectivity for fuel
  xt::xtensor<int, 4> clad_cells = vpin.clad_connectivity();
  xt::xtensor<int, 4> clad_cell_types = xt::view(clad_cells, xt::all(), xt::all(), xt::all(), xt::range(0,1));
  clad_cells = xt::view(clad_cells, xt::all(), xt::all(), xt::all(), xt::range(1, _));
  // adjust cladding connctivity by the number of existing fuel points
  xt::view(clad_cells, xt::all(), xt::all(), xt::all(), xt::range(1,_)) += fuel_points;


  // open vtk file
  std::ofstream fh(filename, std::ofstream::out);

  fh << "# vtk DataFile Version 2.0\n";
  fh << "No comment\nASCII\nDATASET UNSTRUCTURED_GRID\n";
  fh << "POINTS " << total_points << " float\n";
  xt::xtensor<double, 1> points_flat = xt::flatten(pin_points);
  for (auto p = points_flat.begin(); p != points_flat.end(); p+=3) {
    fh << *p << " " << *(p+1) << " " << *(p+2) << "\n";
  }

  points_flat = xt::flatten(clad_points);
  for (auto p = points_flat.begin(); p != points_flat.end(); p+=3) {
    fh << *p << " " << *(p+1) << " " << *(p+2) << "\n";
  }

  // write connectivity header
  fh << "\nCELLS " << num_mesh_elements << " " << num_entries << "\n";

  // write mesh element connectivity
  xt::xtensor<int, 1> cells_flat = xt::flatten(cells);
  int conn_size = HEX_SIZE_ + 1;
  for (auto c = cells_flat.begin(); c != cells_flat.end(); c += conn_size) {
    for (int i = 0; i < conn_size; i++) {
      auto val = *(c+i);
      if (val >= 0) { fh << val << " "; }
    }
    fh << "\n";
  }

  // write cladding connectivity
  cells_flat = xt::flatten(clad_cells);
  conn_size = HEX_SIZE_ + 1;
  for (auto c = cells_flat.begin(); c != cells_flat.end(); c += conn_size) {
    for (int i = 0; i < conn_size; i++) {
      auto val = *(c+i);
      if (val >= 0) { fh << val << " "; }
    }
    fh << "\n";
  }

  // write cell types
  fh << "\nCELL_TYPES " << num_mesh_elements << "\n";
  for (auto v : cell_types) {
    fh << v << "\n";
  }

  // write cell types
  for (auto v : clad_cell_types) {
    fh << v << "\n";
  }

  // data header
  fh << "CELL_DATA " << num_mesh_elements << "\n";

  // number of radial rings in a pin, determines
  // how many times to write a data value
  int num_rings;
  // radial rings in the fuel
  num_rings  = (r_grid_fuel_.size() - 1);
  // radial rings in the cladding
  num_rings += (r_grid_clad_.size() - 1);
  std::cout << "NUM RINGS: " << num_rings << std::endl;

  // temperature data
  fh << "SCALARS TEMPERATURE double 1\n";
  fh << "LOOKUP_TABLE default\n";
  for (int i = 0; i < z_.size() - 1; i++) {
    for (int j = 0; j < num_rings; j++) {
      for (int k = 0; k < radial_resolution; k++) {
        fh << temperature_(0, i, j) << "\n";
      }
    }
  }

  // fission source data
  fh << "SCALARS SOURCE double 1\n";
  fh << "LOOKUP_TABLE default\n";
  for (int i = 0; i < z_.size() - 1; i++) {
    for (int j = 0; j < num_rings; j++) {
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
                                   xt::xtensor<double, 1> c_grid,
                                   int t_res)
  : x_(x), y_(y), z_grid_(z_grid), r_grid_(r_grid), c_grid_(c_grid), t_res_(t_res)
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

xt::xtensor<double, 2> SurrogateHeatDriver::VisualizationPin::create_ring(double radius,
                                                                          int t_resolution) {
  xt::xtensor<double, 1> x = xt::zeros<double>({t_resolution});
  xt::xtensor<double, 1> y = xt::zeros<double>({t_resolution});

  xt::xtensor<double, 1> theta;
  theta = xt::linspace<double>(0., 2.*openmc::PI, t_resolution);

  x = radius * xt::cos(theta);
  y = radius * xt::sin(theta);

  xt::xtensor<double, 2> out = xt::zeros<double>({2, t_resolution});

  xt::view(out, 0, xt::all()) = x;
  xt::view(out, 1, xt::all()) = y;

  return out;
}

xt::xtensor<double, 3> SurrogateHeatDriver::VisualizationPin::fuel_points() {

  xt::xarray<double> pnts_out = xt::zeros<double>({axial_divs_ + 1,
                                                  points_per_plane_,
                                                  3});

  xt::xarray<double> x = xt::zeros<double>({radial_divs_, t_res_});
  xt::xarray<double> y = xt::zeros<double>({radial_divs_, t_res_});

  xt::xtensor<double, 1> theta;
  theta = xt::linspace<double>(0., 2.*openmc::PI, t_res_ + 1);

  // first point is the pin center, start at one
  for(int i = 0; i < radial_divs_; i++) {
    double ring_rad = r_grid_(i);
    xt::xtensor<double, 2> ring = create_ring(ring_rad, t_res_);
    xt::view(x, i, xt::all()) = xt::view(ring, 0, xt::all());
    xt::view(y, i, xt::all()) = xt::view(ring, 1, xt::all());
  }
  // flatten x,y point arrays
  x = xt::flatten(x);
  y = xt::flatten(y);


  for (int i = 0; i < z_grid_.size(); i++) {
    // set all but the center point
    xt::view(pnts_out, i, xt::range(1,_), 0) = x;
    xt::view(pnts_out, i, xt::range(1,_), 1) = y;
    xt::view(pnts_out, i, xt::all(), 2) = z_grid_[i];
  }

  // translate
  xt::view(pnts_out, xt::all(), xt::all(), 0) += x_;
  xt::view(pnts_out, xt::all(), xt::all(), 1) += y_;

  return pnts_out;
}

xt::xtensor<int, 4> SurrogateHeatDriver::VisualizationPin::fuel_connectivity() {
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

xt::xtensor<double, 3> SurrogateHeatDriver::VisualizationPin::clad_points() {

  int clad_points_per_plane = c_grid_.size()*t_res_;
  int clad_divs_ = c_grid_.size();
  xt::xarray<double> pnts_out = xt::zeros<double>({axial_divs_ + 1,
                                                  clad_points_per_plane,
                                                  3});

  xt::xarray<double> x = xt::zeros<double>({clad_divs_, t_res_});
  xt::xarray<double> y = xt::zeros<double>({clad_divs_, t_res_});

  xt::xtensor<double, 1> theta;
  theta = xt::linspace<double>(0., 2.*openmc::PI, t_res_ + 1);

  for(int i = 0; i < clad_divs_; i++) {
    double ring_rad = c_grid_(i);
    xt::xtensor<double, 2> ring = create_ring(ring_rad, t_res_);
    xt::view(x, i, xt::all()) = xt::view(ring, 0, xt::all());
    xt::view(y, i, xt::all()) = xt::view(ring, 1, xt::all());
  }
  // flatten x,y point arrays
  x = xt::flatten(x);
  y = xt::flatten(y);


  for (int i = 0; i < z_grid_.size(); i++) {
    // set all but the center point
    xt::view(pnts_out, i, xt::all(), 0) = x;
    xt::view(pnts_out, i, xt::all(), 1) = y;
    xt::view(pnts_out, i, xt::all(), 2) = z_grid_[i];
  }

  // translate
  xt::view(pnts_out, 0) += x_;
  xt::view(pnts_out, 1) += y_;

  return pnts_out;
}

xt::xtensor<int, 4> SurrogateHeatDriver::VisualizationPin::clad_connectivity() {
  int clad_points_per_plane = c_grid_.size()*t_res_;
  int clad_divs_ = c_grid_.size() - 1;

  // size output array
  xt::xtensor<int, 4> cells_out = xt::zeros<int>({axial_divs_,
                                                  clad_divs_,
                                                  t_res_,
                                                  HEX_SIZE_ + 2});

  // generate a base layer to be extended in Z
  xt::xtensor<int, 3> base = xt::zeros<int>({clad_divs_, t_res_, HEX_SIZE_});

  /// OUTER RINGS \\\

  xt::xtensor<int, 2> radial_base = xt::zeros<int>({t_res_, HEX_SIZE_});

  // setup connectivity of the first layer
  xt::view(radial_base, xt::all(), 0) = xt::arange(0, t_res_);
  xt::view(radial_base, xt::all(), 1) = xt::arange(1, t_res_ + 1);
  xt::view(radial_base, xt::all(), 2) = xt::view(radial_base, xt::all(), 1) + t_res_;
  xt::view(radial_base, xt::all(), 3) = xt::view(radial_base, xt::all(), 0) + t_res_;
  xt::view(radial_base, t_res_ - 1, 1) = 1;
  xt::view(radial_base, t_res_ - 1, 2) = t_res_ + 1;

  // copy connectivity of the first layer
  xt::view(radial_base, xt::all(), xt::range(4,8)) +=
    xt::view(radial_base, xt::all(), xt::range(0,4));
  // shift connectivity down one layer
  xt::view(radial_base, xt::all(), xt::range(4,8)) += clad_points_per_plane;

  // other rings
  xt::strided_view(base, {xt::range(0, _), xt::ellipsis()}) = radial_base;
  for (int i = 0; i < clad_divs_; i++) {
    // extend based on the starting index of the ring
    int start_idx = i * t_res_;
    xt::view(base, i, xt::all(), xt::all()) += start_idx;
  }

  // the rest are hexes
  xt::view(cells_out, xt::all(), xt::all(), xt::all(), 0) = HEX_TYPE_;
  xt::view(cells_out, xt::all(), xt::all(), xt::all(), 1) = HEX_SIZE_;

  // set all axial divs using base
  for(int i = 0; i < axial_divs_; i++) {
    // set layer and increment connectivity by number of points in axial div
    xt::view(cells_out, i, xt::all(), xt::all(), xt::range(2, 10)) = base;
    base += clad_points_per_plane;
  }

  return cells_out;
}


} // namespace stream
