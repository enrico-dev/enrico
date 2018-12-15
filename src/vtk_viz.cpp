

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

SurrogateToVtk::SurrogateToVtk(const SurrogateHeatDriver* surrogate_ptr) :
sgate(surrogate_ptr) {

  // create a pin
   radial_res = 5;

   n_axial_sections = sgate->z_.size() - 1;
   n_axial_points = sgate->z_.size();
   n_radial_fuel_sections = sgate->r_grid_fuel_.size() - 1;
   n_radial_clad_sections = sgate->r_grid_clad_.size() - 1;
   n_radial_sections = n_radial_fuel_sections + n_radial_clad_sections;
   n_sections_per_plane = n_radial_sections * radial_res;
  // Calculate some necessary values ahead of time
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

void SurrogateToVtk::write_vtk() {

  std::ofstream fh("magnolia.vtk", std::ofstream::out);

  fh << "# vtk DataFile Version 2.0\n";
  fh << "No comment\nASCII\nDATASET UNSTRUCTURED_GRID\n";
  fh << "POINTS " << n_fuel_points << " float\n";

  xt::xtensor<double, 1> pnts = points();
  for (auto val = pnts.begin(); val != pnts.end(); val+=3) {
    fh << *val << " " << *(val+1) << " " << *(val+2) << "\n";
  }
  // write connectivity header
  fh << "\nCELLS " << n_fuel_elements << " " << n_fuel_entries_per_plane * n_axial_sections << "\n";

  xt::xtensor<int, 1> connectivity = conn();
  for (auto val = connectivity.begin(); val != connectivity.end(); val += CONN_STRIDE_) {
    for (int i = 0; i < CONN_STRIDE_; i++) {
      auto v = *(val+i);
      if (v >=0) { fh << v << " "; }
    }
    fh << "\n";
  }

  fh << "\nCELL_TYPES " << n_fuel_elements << "\n";
  xt::xtensor<int ,1> cell_types = types();
  for (auto v : cell_types) {
      fh << v << "\n";
  }

}

xt::xtensor<double, 1> SurrogateToVtk::points() {
  xt::xtensor<double, 3> fuel_pnts = fuel_points();
  xt::xtensor<double, 3> clad_pnts = clad_points();
  xt::xtensor<double, 1> points = xt::flatten(fuel_pnts);
  return points;
}

xt::xtensor<double, 3> SurrogateToVtk::fuel_points() {

  xt::xarray<double> pnts_out = xt::zeros<double>({n_axial_points,
                                                  fuel_points_per_plane,
                                                  3});

  xt::xarray<double> x = xt::zeros<double>({n_radial_fuel_sections, radial_res});
  xt::xarray<double> y = xt::zeros<double>({n_radial_fuel_sections, radial_res});

  for(int i = 0; i < n_radial_fuel_sections; i++) {
    // first radius is zero, start at one
    double ring_rad = sgate->r_grid_fuel_(i + 1);
    xt::xtensor<double, 2> ring = create_ring(ring_rad, radial_res);
    xt::view(x, i, xt::all()) = xt::view(ring, 0, xt::all());
    xt::view(y, i, xt::all()) = xt::view(ring, 1, xt::all());
  }
  // flatten x,y point arrays
  x = xt::flatten(x);
  y = xt::flatten(y);

  for (int i = 0; i < n_axial_points; i++) {
    // set all but the center point
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

  xt::xarray<double> pnts_out = xt::zeros<double>({n_axial_points,
                                                  clad_points_per_plane,
                                                  3});

  xt::xarray<double> x = xt::zeros<double>({n_radial_clad_sections + 1, radial_res});
  xt::xarray<double> y = xt::zeros<double>({n_radial_clad_sections + 1, radial_res});

  for(int i = 0; i < n_radial_clad_sections; i++) {
    double ring_rad = sgate->r_grid_clad_(i);
    xt::xtensor<double, 2> ring = create_ring(ring_rad, radial_res);
    xt::view(x, i, xt::all()) = xt::view(ring, 0, xt::all());
    xt::view(y, i, xt::all()) = xt::view(ring, 1, xt::all());
  }
  // flatten x,y point arrays
  x = xt::flatten(x);
  y = xt::flatten(y);

  for (int i = 0; i < n_axial_points; i++) {
    // set all but the center point
    xt::view(pnts_out, i, xt::all(), 0) = x;
    xt::view(pnts_out, i, xt::all(), 1) = y;
    xt::view(pnts_out, i, xt::all(), 2) = sgate->z_[i];
  }

  // translate
  xt::view(pnts_out, 0) += sgate->pin_centers_(0,0);
  xt::view(pnts_out, 1) += sgate->pin_centers_(0,1);

  return pnts_out;
}

xt::xtensor<int, 1> SurrogateToVtk::conn() {
  xt::xtensor<int, 4> f_conn = fuel_conn();
  xt::xtensor<int , 1> out = xt::flatten(f_conn);
  return out;
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


  /// INNER RING \\\

  xt::xtensor<int, 2> inner_base = xt::zeros<int>({radial_res, HEX_SIZE_});

  // cell connectivity for the first z level
  xt::view(inner_base, xt::all(), 1) = xt::arange(1, radial_res + 1);
  xt::view(inner_base, xt::all(), 2) = xt::arange(2, radial_res + 2);
  // adjust last cell for perioic condition
   xt::view(inner_base, radial_res - 1, 2) = 1;

  // copy connectivity of first layer to the second
  xt::view(inner_base, xt::all(), xt::range(3,6)) =
    xt::view(inner_base, xt::all(), xt::range(0,3));
  // shift connectivity down one layer
  xt::strided_view(inner_base, {xt::all(), xt::range(3,6)}) += fuel_points_per_plane;

  // set the inner_base
   xt::view(base, 0, xt::all(), xt::all()) = inner_base;

  /// OUTER RINGS \\\

  xt::xtensor<int, 2> radial_base = xt::zeros<int>({radial_res, HEX_SIZE_});

  // setup connectivity of the first layer
  xt::view(radial_base, xt::all(), 0) = xt::arange(1, radial_res + 1);
  xt::view(radial_base, xt::all(), 1) = xt::arange(2, radial_res + 2);
  xt::view(radial_base, xt::all(), 2) = xt::view(radial_base, xt::all(), 1) + radial_res;
  xt::view(radial_base, xt::all(), 3) = xt::view(radial_base, xt::all(), 0) + radial_res;
  xt::view(radial_base, radial_res - 1, 1) = 1;
  xt::view(radial_base, radial_res - 1, 2) = radial_res + 1;

  // copy connectivity of the first layer
  xt::view(radial_base, xt::all(), xt::range(4,8)) +=
    xt::view(radial_base, xt::all(), xt::range(0,4));
  // shift connectivity down one layer
  xt::view(radial_base, xt::all(), xt::range(4,8)) += fuel_points_per_plane;

  // other rings
  xt::strided_view(base, {xt::range(1, _), xt::ellipsis()}) = radial_base;
  for (int i = 1; i < n_radial_fuel_sections; i++) {
    // extend based on the starting index of the ring
    int start_idx = (i-1) * radial_res;
    xt::view(base, i, xt::all(), xt::all()) += start_idx;
  }

  // set all axial divs using base
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

xt::xtensor<int, 1> SurrogateToVtk::types() {
  xt::xtensor<int, 3> ftypes = fuel_types();
  return xt::flatten(ftypes);
}

xt::xtensor<int, 3> SurrogateToVtk::fuel_types() {

  xt::xtensor<int, 3> types_out = xt::zeros<int>({n_axial_sections,
                                                  n_radial_fuel_sections,
                                                  radial_res});

  // the inner ring is wedges
  xt::view(types_out, xt::all(), 0, xt::all()) = WEDGE_TYPE_;

  // the rest are hexes
  xt::view(types_out, xt::all(), xt::range(1, _), xt::all()) = HEX_TYPE_;

  return types_out;
}


} // stream
