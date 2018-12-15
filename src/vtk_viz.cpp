

#include "stream/vtk_viz.h"
#include "xtensor/xtensor.hpp"

#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xview.hpp"

#include "openmc/constants.h"

// some constant values
const int WEDGE_TYPE_ = 13;
const int WEDGE_SIZE_ = 6;
const int HEX_TYPE_ = 12;
const int HEX_SIZE_ = 8;
const int INVALID_CONN_ = -1;
const int CONN_STRIDE_ = HEX_SIZE_ + 1;

namespace stream {

using namespace xt::placeholders;

xt::xtensor<double, 2> create_ring(double radius,
                                   int t_resolution) {
  xt::xtensor<double, 1> x = xt::zeros<double>({t_resolution});
  xt::xtensor<double, 1> y = xt::zeros<double>({t_resolution});

  xt::xtensor<double, 1> theta;
  theta = xt::linspace<double>(0., 2.*openmc::PI, t_resolution + 1);
  theta = xt::view(theta, xt::range(0, -1)); // remove last value
  x = radius * xt::cos(theta);
  y = radius * xt::sin(theta);

  xt::xtensor<double, 2> out = xt::zeros<double>({2, t_resolution});

  xt::view(out, 0, xt::all()) = x;
  xt::view(out, 1, xt::all()) = y;

  return out;
}

xt::xtensor<int, 2> hex_ring(int start_idx, int resolution, int z_shift) {
  xt::xtensor<int, 2> out = xt::zeros<int>({resolution, HEX_SIZE_});

  // set connectivity of the first z-layer
  xt::view(out, xt::all(), 0) = xt::arange(start_idx, resolution + start_idx);
  xt::view(out, xt::all(), 1) = xt::arange(start_idx + 1, resolution + start_idx + 1);
  xt::view(out, xt::all(), 2) = xt::view(out, xt::all(), 1) + resolution;
  xt::view(out, xt::all(), 3) = xt::view(out, xt::all(), 0) + resolution;
  // adjust last point id for periodic condition on both layers, using hexes now
  xt::view(out, resolution - 1, 1) = 1;
  xt::view(out, resolution - 1, 2) = resolution + 1;

  // copy connectivity of the first layer to the second
  // and shift by number of points
  xt::view(out, xt::all(), xt::range(4,HEX_SIZE_)) +=
    xt::view(out, xt::all(), xt::range(0,4));
  // shift connectivity down one layer
  xt::view(out, xt::all(), xt::range(4,HEX_SIZE_)) += z_shift;

  return out;
}

SurrogateToVtk::SurrogateToVtk(const SurrogateHeatDriver* surrogate_ptr, int t_res) :
sgate(surrogate_ptr), radial_res(t_res) {
  // Set some necessary values ahead of time
   n_axial_sections = sgate->z_.size() - 1;
   n_axial_points = sgate->z_.size();
   n_radial_fuel_sections = sgate->r_grid_fuel_.size() - 1;
   n_radial_clad_sections = sgate->r_grid_clad_.size() - 1;
   n_radial_sections = n_radial_fuel_sections + n_radial_clad_sections;
   n_sections_per_plane = n_radial_sections * radial_res;
  // fuel points
  fuel_points_per_plane  = n_radial_fuel_sections * radial_res + 1;
  n_fuel_points = fuel_points_per_plane * sgate->z_.size();
  // cladding points
  clad_points_per_plane = sgate->r_grid_clad_.size() * radial_res;
  n_clad_points = clad_points_per_plane * sgate->r_grid_clad_.size();

  points_per_plane = fuel_points_per_plane + clad_points_per_plane;
  // total number of points in pin
  total_points = points_per_plane * sgate->z_.size();

  // fuel elements, entries
  n_mesh_elements = n_sections_per_plane * n_axial_sections;
  n_fuel_elements = n_radial_fuel_sections * radial_res * n_axial_sections;
  n_clad_elements = n_radial_clad_sections * radial_res * n_axial_sections;
  // wedge regions
  n_fuel_entries_per_plane = radial_res * (WEDGE_SIZE_ + 1);
  // other radial regions
  n_fuel_entries_per_plane += radial_res * (HEX_SIZE_ + 1) * (n_radial_fuel_sections - 1);
  n_clad_entries_per_plane = radial_res * (HEX_SIZE_ + 1) * n_radial_clad_sections;
  n_entries_per_plane = n_fuel_entries_per_plane + n_clad_entries_per_plane;
  // cladding elements
  n_entries = n_entries_per_plane * n_axial_sections;
}

void SurrogateToVtk::write_vtk(std::string filename) {

  std::ofstream fh(filename, std::ofstream::out);

  // write vtk header
  fh << "# vtk DataFile Version 2.0\n";
  fh << "No comment\nASCII\nDATASET UNSTRUCTURED_GRID\n";

  /// POINTS \\\

  fh << "POINTS " << total_points << " float\n";
  xt::xtensor<double, 1> pnts = points();
  for (auto val = pnts.begin(); val != pnts.end(); val+=3) {
    fh << *val << " " << *(val+1) << " " << *(val+2) << "\n";
  }

  /// ELEMENT CONNECTIVITY \\\

  fh << "\nCELLS " << n_mesh_elements << " " << n_entries << "\n";
  xt::xtensor<int, 1> connectivity = conn();
  for (auto val = connectivity.begin(); val != connectivity.end(); val += CONN_STRIDE_) {
    for (int i = 0; i < CONN_STRIDE_; i++) {
      auto v = *(val+i);
      if (v >=0) { fh << v << " "; }
    }
    fh << "\n";
  }

  /// MESH ELEMENT TYPES \\\

  fh << "\nCELL_TYPES " << n_mesh_elements << "\n";
  xt::xtensor<int ,1> cell_types = types();
  for (auto v : cell_types) {
      fh << v << "\n";
  }

  /// WRITE DATA \\\
  // fuel mesh elements are written first, followed by cladding elements
  // the data needs to be written in a similar matter

  // for each radial section, the data point for that radial ring
  // is repeated radial_res times
  fh << "CELL_DATA " << n_mesh_elements << "\n";

  // temperature data
  fh << "SCALARS TEMPERATURE double 1\n";
  fh << "LOOKUP_TABLE default\n";
  // write all fuel data first
  for (int i = 0; i < n_axial_sections; i++) {
    for (int j = 0; j < n_radial_fuel_sections; j++) {
      for (int k = 0; k < radial_res; k++) {
        fh << sgate->temperature_(0, i, j) << "\n";
      }
    }
  }
  // then write cladding data
  for (int i = 0; i < n_axial_sections; i++) {
    for (int j = 0; j < n_radial_clad_sections; j++) {
      for (int k = 0; k < radial_res; k++) {
        fh << sgate->temperature_(0, i, j + n_radial_fuel_sections) << "\n";
      }
    }
  }

  // source data
  fh << "SCALARS SOURCE double 1\n";
  fh << "LOOKUP_TABLE default\n";
  // write all fuel data first
  for (int i = 0; i < n_axial_sections; i++) {
    for (int j = 0; j < n_radial_fuel_sections; j++) {
      for (int k = 0; k < radial_res; k++) {
        fh << sgate->source_(0, i, j) << "\n";
      }
    }
  }
  // then write the cladding data
  for (int i = 0; i < n_axial_sections; i++) {
    for (int j = 0; j < n_radial_clad_sections; j++) {
      for (int k = 0; k < radial_res; k++) {
        fh << sgate->source_(0, i, j + n_radial_fuel_sections) << "\n";
      }
    }
  }

  // close the file
  fh.close();

} // write vtk

xt::xtensor<double, 1> SurrogateToVtk::points() {
  // get fuel points
  xt::xtensor<double, 3> fuel_pnts = fuel_points();
  // get cladding points
  xt::xtensor<double, 3> clad_pnts = clad_points();
  // concatenate in 1-D and return
  xt::xtensor<double, 1> points = xt::concatenate(xt::xtuple(
    xt::flatten(fuel_pnts), xt::flatten(clad_pnts)));
  return points;
}

xt::xtensor<double, 3> SurrogateToVtk::fuel_points() {
  // array to hold all point data
  xt::xarray<double> pnts_out = xt::zeros<double>({n_axial_points,
                                                  fuel_points_per_plane,
                                                  3});

  // x and y for a single axial plane in the rod
  xt::xarray<double> x = xt::zeros<double>({n_radial_fuel_sections, radial_res});
  xt::xarray<double> y = xt::zeros<double>({n_radial_fuel_sections, radial_res});

  // generate x/y points for each raidal section
  for(int i = 0; i < n_radial_fuel_sections; i++) {
    // first radius in r_grid_fuel_ is zero, start at one
    double ring_rad = sgate->r_grid_fuel_(i + 1);
    // create a unit ring scaled by radius
    xt::xtensor<double, 2> ring = create_ring(ring_rad, radial_res);
    // set x/y values for this ring
    xt::view(x, i, xt::all()) = xt::view(ring, 0, xt::all());
    xt::view(y, i, xt::all()) = xt::view(ring, 1, xt::all());
  }

  // flatten radial point arrays
  x = xt::flatten(x);
  y = xt::flatten(y);

  // set the point values for each axial plane
  for (int i = 0; i < n_axial_points; i++) {
    // set all but the center point, which is always (0,0,z)
    xt::view(pnts_out, i, xt::range(1,_), 0) = x;
    xt::view(pnts_out, i, xt::range(1,_), 1) = y;
    xt::view(pnts_out, i, xt::all(), 2) = sgate->z_(i);
  }

  // translate
  xt::view(pnts_out, xt::all(), xt::all(), 0) += sgate->pin_centers_(0,0);
  xt::view(pnts_out, xt::all(), xt::all(), 1) += sgate->pin_centers_(0,1);

  return pnts_out;
}

xt::xtensor<double, 3> SurrogateToVtk::clad_points() {
  // array to hold all point data
  xt::xarray<double> pnts_out = xt::zeros<double>({n_axial_points,
                                                  clad_points_per_plane,
                                                  3});
  // x and y for a single axial plane in the rod
  xt::xarray<double> x = xt::zeros<double>({n_radial_clad_sections + 1, radial_res});
  xt::xarray<double> y = xt::zeros<double>({n_radial_clad_sections + 1, radial_res});

  // generate x/y points for each radial ring
  for(int i = 0; i < n_radial_clad_sections + 1; i++) {
    double ring_rad = sgate->r_grid_clad_(i);
    // create a unit ring scaled by radius
    xt::xtensor<double, 2> ring = create_ring(ring_rad, radial_res);
    // set x/y values for this ring
    xt::view(x, i, xt::all()) = xt::view(ring, 0, xt::all());
    xt::view(y, i, xt::all()) = xt::view(ring, 1, xt::all());
  }

  // flatten x,y point arrays
  x = xt::flatten(x);
  y = xt::flatten(y);

  for (int i = 0; i < n_axial_points; i++) {
    // set all points, no center point for cladding
    xt::view(pnts_out, i, xt::all(), 0) = x;
    xt::view(pnts_out, i, xt::all(), 1) = y;
    xt::view(pnts_out, i, xt::all(), 2) = sgate->z_[i];
  }

  // translate to pin center
  xt::view(pnts_out, 0) += sgate->pin_centers_(0,0);
  xt::view(pnts_out, 1) += sgate->pin_centers_(0,1);

  return pnts_out;
}

xt::xtensor<int, 1> SurrogateToVtk::conn() {
  // get the fuel connectivity
  xt::xtensor<int, 4> f_conn = fuel_conn();
  // get the cladding connectivity
  xt::xtensor<int, 4> c_conn = clad_conn();
  // adjust the cladding connectivity by
  // the number of points in the fuel mesh
  xt::view(c_conn, xt::all(), xt::all(), xt::all(), xt::range(1, _)) += fuel_points_per_plane * n_axial_points;
  // concatenate and return 1-D form
  xt::xtensor<int, 1> conn_out = xt::concatenate(xt::xtuple(
    xt::flatten(f_conn), xt::flatten(c_conn)));
  return conn_out;
}

xt::xtensor<int, 4> SurrogateToVtk::fuel_conn() {
  // size output array
  xt::xtensor<int, 4> cells_out = xt::zeros<int>({n_axial_sections,
                                                  n_radial_fuel_sections,
                                                  radial_res,
                                                  HEX_SIZE_ + 1});

  // generate a base layer to be extended in Z
  xt::xtensor<int, 3> base = xt::zeros<int>({n_radial_fuel_sections,
                                             radial_res,
                                             HEX_SIZE_});


  /// INNERMOST RING (WEDGES) \\\

  xt::xtensor<int, 2> inner_base = xt::zeros<int>({radial_res, HEX_SIZE_});

  xt::view(inner_base, xt::all(), 1) = xt::arange(1, radial_res + 1);
  xt::view(inner_base, xt::all(), 2) = xt::arange(2, radial_res + 2);
  // adjust last point id for perioic condition
   xt::view(inner_base, radial_res - 1, 2) = 1;

  // copy connectivity of first z-layer to the second and shift
  // by the number of points in a plane
  xt::view(inner_base, xt::all(), xt::range(3,6)) =
    xt::view(inner_base, xt::all(), xt::range(0,3));
  xt::strided_view(inner_base, {xt::all(), xt::range(3,6)}) += fuel_points_per_plane;

  // set values for the innermost ring
  xt::view(base, 0, xt::all(), xt::all()) = inner_base;

  /// OUTER RINGS (HEXES) \\\

  // create a radial base starting at point zero
  // with a shift between the first and second layer of connectivity
  // equal to the number of fuel points in a plane
  xt::xtensor<int, 2> radial_base = hex_ring(1, radial_res, fuel_points_per_plane);

  // set the other rings by shifting the initial base
  // by radial_res for each ring
  for (int i = 1; i < n_radial_fuel_sections; i++) {
    xt::view(base, i, xt::all(), xt::all()) = radial_base;
    radial_base += radial_res;
  }

  // set all axial divs using the base connectivity
  // and shifting by the number of points in a plane
  for (int j = 0; j < n_axial_sections; j++) {
    // set layer and increment connectivity by number of points in axial div
    xt::view(cells_out, j, xt::all(), xt::all(), xt::range(1, _)) = base;
    base += fuel_points_per_plane;
  }

  // innermost ring is always wedges
  xt::view(cells_out, xt::all(), 0, xt::all(), 0) = WEDGE_SIZE_;
  // the reset are hexes
  xt::view(cells_out, xt::all(), xt::range(1, _), xt::all(), 0) = HEX_SIZE_;
  // first ring should be wedges only, invalidate last two entries
  xt::view(cells_out, xt::all(), 0, xt::all(), xt::range(7,_)) = INVALID_CONN_;

  return cells_out;
}

xt::xtensor<int, 4> SurrogateToVtk::clad_conn() {

  // size output array
  xt::xtensor<int, 4> cells_out = xt::zeros<int>({n_axial_sections,
                                                  n_radial_clad_sections,
                                                  radial_res,
                                                  HEX_SIZE_ + 1});


  // base layer to be extended in Z
  xt::xtensor<int, 3> base = xt::zeros<int>({n_radial_clad_sections,
                                             radial_res,
                                             HEX_SIZE_});


  /// ELEMENT RINGS (HEXES) \\\

  // create a radial base starting at point zero
  // with a shift between the first and second layer of connectivity
  // equal to the number of cladding points in a plane
  xt::xtensor<int, 2> radial_base = hex_ring(0, radial_res, clad_points_per_plane);

  // set the other rings by shifting the initial base
  // by radial_res for each ring
  for (int i = 0; i < n_radial_clad_sections; i++) {
    xt::view(base, i, xt::all(), xt::all()) = radial_base;
    radial_base += radial_res;
  }

  // set all axial divs using base for the first layer and
  // shift by the number of points in a plane
  for(int i = 0; i < n_axial_sections; i++) {
    // set layer and increment connectivity by number of points in axial div
    xt::view(cells_out, i, xt::all(), xt::all(), xt::range(1, _)) = base;
    base += clad_points_per_plane;
  }

  // all are hexes
  xt::view(cells_out, xt::all(), xt::all(), xt::all(), 0) = HEX_SIZE_;

  return cells_out;
}

xt::xtensor<int, 1> SurrogateToVtk::types() {
  // get fuel types
  xt::xtensor<int, 3> ftypes = fuel_types();
  // get the cladding types
  xt::xtensor<int, 3> ctypes = clad_types();
  // concatenate and return 1-D form
  xt::xtensor<int, 1> types_out = xt::concatenate(xt::xtuple(
    xt::flatten(ftypes), xt::flatten(ctypes)));
  return types_out;
}

xt::xtensor<int, 3> SurrogateToVtk::fuel_types() {
  // size the output array
  xt::xtensor<int, 3> types_out = xt::zeros<int>({n_axial_sections,
                                                  n_radial_fuel_sections,
                                                  radial_res});

  // the inner ring is wedges
  xt::view(types_out, xt::all(), 0, xt::all()) = WEDGE_TYPE_;
  // the rest are hexes
  xt::view(types_out, xt::all(), xt::range(1, _), xt::all()) = HEX_TYPE_;

  return types_out;
}

xt::xtensor<int, 3> SurrogateToVtk::clad_types() {
  // size output array
  xt::xtensor<int, 3> clad_types_out = xt::zeros<int>({n_axial_sections,
                                                      n_radial_clad_sections,
                                                      radial_res});
  // all elements are hexes
  clad_types_out = xt::full_like(clad_types_out, HEX_TYPE_);

  return clad_types_out;
}

} // stream
