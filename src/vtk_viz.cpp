#include <cmath>

#include "enrico/vtk_viz.h"

#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "openmc/constants.h"

// some constant values
const int WEDGE_TYPE_ = 13;
const size_t WEDGE_SIZE_ = 6;
const int HEX_TYPE_ = 12;
const size_t HEX_SIZE_ = 8;
const int INVALID_CONN_ = -1;
const size_t CONN_STRIDE_ = HEX_SIZE_ + 1;

namespace enrico {

using std::ofstream;
using xt::xtensor;
using xt::placeholders::_;

xtensor<double, 2> create_ring(double radius, size_t t_resolution)
{
  xtensor<double, 1> x({t_resolution}, 0.0);
  xtensor<double, 1> y({t_resolution}, 0.0);

  xtensor<double, 1> theta;
  theta = xt::linspace<double>(0., 2. * M_PI, t_resolution + 1);
  theta = xt::view(theta, xt::range(0, -1)); // remove last value

  x = radius * xt::cos(theta);
  y = radius * xt::sin(theta);

  xtensor<double, 2> out({2, (size_t)t_resolution}, 0.0);

  xt::view(out, 0, xt::all()) = x;
  xt::view(out, 1, xt::all()) = y;

  return out;
}

xtensor<int, 2> hex_ring(size_t start_idx, size_t resolution, size_t z_shift)
{
  xtensor<int, 2> out({resolution, HEX_SIZE_}, 0);

  // set connectivity of the first z-layer
  // first two points - along the inner radial ring
  xt::view(out, xt::all(), 0) = xt::arange(start_idx, resolution + start_idx);
  xt::view(out, xt::all(), 1) = xt::arange(start_idx + 1, resolution + start_idx + 1);

  // second two poitns - along the outer radial ring, going back toward the starting point
  xt::view(out, xt::all(), 2) = xt::view(out, xt::all(), 1) + resolution;
  xt::view(out, xt::all(), 3) = xt::view(out, xt::all(), 0) + resolution;

  // adjust last point id for periodic condition on both layers, using hexes now
  xt::view(out, resolution - 1, 1) = start_idx;
  xt::view(out, resolution - 1, 2) = resolution + start_idx;

  // copy connectivity of the first layer to the second
  // and shift by number of points
  xt::view(out, xt::all(), xt::range(4, HEX_SIZE_)) +=
    xt::view(out, xt::all(), xt::range(0, 4));

  // shift connectivity down one layer
  xt::view(out, xt::all(), xt::range(4, HEX_SIZE_)) += z_shift;

  return out;
}

SurrogateVtkWriter::SurrogateVtkWriter(const SurrogateHeatDriver& surrogate_ref,
                                       size_t t_res,
                                       const std::string& regions_to_write,
                                       const std::string& data_to_write)
  : surrogate_(surrogate_ref)
  , azimuthal_res_(t_res)
{

  // read data specs
  data_out_ = VizDataType::all;
  if ("all" == data_to_write) {
    data_out_ = VizDataType::all;
  } else if ("source" == data_to_write) {
    data_out_ = VizDataType::source;
  } else if ("temp" == data_to_write || "temperature" == data_to_write) {
    data_out_ = VizDataType::temp;
  } else if ("density" == data_to_write) {
    data_out_ = VizDataType::density;
  } else {
    // invalid user input
    Expects(false);
  }

  output_includes_temp_ =
    (data_out_ == VizDataType::all) || (data_out_ == VizDataType::temp);

  output_includes_density_ =
    (data_out_ == VizDataType::all) || (data_out_ == VizDataType::density);

  output_includes_source_ =
    (data_out_ == VizDataType::all) || (data_out_ == VizDataType::source);

  // read data specs
  regions_out_ = VizRegionType::all;
  if ("all" == regions_to_write) {
    regions_out_ = VizRegionType::all;
  } else if ("solid" == regions_to_write) {
    regions_out_ = VizRegionType::solid;
  } else if ("fluid" == regions_to_write) {
    regions_out_ = VizRegionType::fluid;
  } else {
    // invalid user input
    Expects(false);
  }

  output_includes_fluid_ =
    (regions_out_ == VizRegionType::all) || (regions_out_ == VizRegionType::fluid);

  output_includes_solid_ =
    (regions_out_ == VizRegionType::all) || (regions_out_ == VizRegionType::solid);

  // if the output includes the fluid phase, for simplicity of constructing
  // the wedges, we require the azimuthal resolution to be divisible by the
  // number of channels around the rod
  if (output_includes_fluid_) {
    Expects(azimuthal_res_ % chans_per_rod_ == 0);
  }

  set_number_of_sections();

  set_number_of_points();

  set_number_of_entries();

  // generate a representative set of points and connectivity
  points_ = points();
  conn_ = conn();
  types_ = types();
}

void SurrogateVtkWriter::set_number_of_sections()
{
  n_axial_sections_ = surrogate_.z_.size() - 1;
  n_radial_fuel_sections_ = surrogate_.n_fuel_rings();
  n_radial_clad_sections_ = surrogate_.n_clad_rings();
  n_fluid_sections_ = azimuthal_res_ + 2 * chans_per_rod_;

  if (regions_out_ == VizRegionType::solid) {
    n_sections_ = (n_radial_fuel_sections_ + n_radial_clad_sections_) * azimuthal_res_ *
                  n_axial_sections_;
  } else if (regions_out_ == VizRegionType::fluid) {
    n_sections_ = n_fluid_sections_ * n_axial_sections_;
  } else if (regions_out_ == VizRegionType::all) {
    n_sections_ = (n_radial_fuel_sections_ + n_radial_clad_sections_) * azimuthal_res_ *
                    n_axial_sections_ +
                  n_fluid_sections_ * n_axial_sections_;
  }
}

void SurrogateVtkWriter::set_number_of_points()
{
  n_axial_points_ = surrogate_.z_.size();
  fuel_points_per_plane_ = n_radial_fuel_sections_ * azimuthal_res_ + 1;
  clad_points_per_plane_ = (n_radial_clad_sections_ + 1) * azimuthal_res_;
  fluid_points_per_plane_ = azimuthal_res_ + 8;

  if (regions_out_ == VizRegionType::solid) {
    n_points_ = (fuel_points_per_plane_ + clad_points_per_plane_) * n_axial_points_;
  } else if (regions_out_ == VizRegionType::fluid) {
    n_points_ = fluid_points_per_plane_ * n_axial_points_;
  } else if (regions_out_ == VizRegionType::all) {
    n_points_ =
      (fuel_points_per_plane_ + clad_points_per_plane_ + fluid_points_per_plane_) *
      n_axial_points_;
  }
}

void SurrogateVtkWriter::set_number_of_entries()
{
  n_fluid_entries_per_plane_ = n_fluid_sections_ * (WEDGE_SIZE_ + 1);

  n_fuel_entries_per_plane_ =
    azimuthal_res_ * (WEDGE_SIZE_ + 1) +
    azimuthal_res_ * (HEX_SIZE_ + 1) * (n_radial_fuel_sections_ - 1);

  n_clad_entries_per_plane_ = azimuthal_res_ * (HEX_SIZE_ + 1) * n_radial_clad_sections_;

  if (regions_out_ == VizRegionType::solid) {
    n_entries_ =
      (n_fuel_entries_per_plane_ + n_clad_entries_per_plane_) * n_axial_sections_;
  } else if (regions_out_ == VizRegionType::fluid) {
    n_entries_ = n_fluid_entries_per_plane_ * n_axial_sections_;
  } else if (regions_out_ == VizRegionType::all) {
    n_entries_ = (n_fuel_entries_per_plane_ + n_clad_entries_per_plane_ +
                  n_fluid_entries_per_plane_) *
                 n_axial_sections_;
  }
}

void SurrogateVtkWriter::write(std::string filename)
{
  // open file
  ofstream fh(filename, std::ofstream::out);

  // write vtk header
  write_header(fh);

  // write vertex locations
  write_points(fh);

  // write wedge/hex element connectivity
  write_element_connectivity(fh);

  // write types for wedge/hex elements
  write_element_types(fh);

  // write specified data to the vtk file
  write_data(fh);

  // close the file
  fh.close();

} // write_vtk

void SurrogateVtkWriter::write_header(ofstream& vtk_file)
{
  vtk_file << "# vtk DataFile Version 2.0\n";
  vtk_file << "No comment\nASCII\nDATASET UNSTRUCTURED_GRID\n";
}

void SurrogateVtkWriter::write_points(ofstream& vtk_file)
{
  vtk_file << "POINTS " << surrogate_.n_pins_ * n_points_ << " float\n";

  for (size_t pin = 0; pin < surrogate_.n_pins_; pin++) {
    // translate pin template to pin center
    xtensor<double, 1> pnts =
      points_for_pin(surrogate_.pin_centers_(pin, 0), surrogate_.pin_centers_(pin, 1));

    for (auto val = pnts.cbegin(); val != pnts.cend(); val += 3) {
      vtk_file << *val << " " << *(val + 1) << " " << *(val + 2) << "\n";
    }
  }
}

void SurrogateVtkWriter::write_element_connectivity(ofstream& vtk_file)
{
  // write number of connectivity entries
  vtk_file << "\nCELLS " << surrogate_.n_pins_ * n_sections_ << " "
           << surrogate_.n_pins_ * n_entries_ << "\n";

  for (size_t pin = 0; pin < surrogate_.n_pins_; pin++) {
    // get the connectivity for a given pin, using an
    // offset to get the connectivity values correct
    xtensor<int, 1> conn = conn_for_pin(pin * n_points_);
    // write the connectivity values to file for this pin
    for (auto val = conn.cbegin(); val != conn.cend(); val += CONN_STRIDE_) {
      for (size_t i = 0; i < CONN_STRIDE_; i++) {
        auto v = *(val + i);
        // mask out any negative connectivity values
        if (v != INVALID_CONN_) {
          vtk_file << v << " ";
        }
      }
      vtk_file << "\n";
    }
  }
} // write_element_connectivity

void SurrogateVtkWriter::write_element_types(ofstream& vtk_file)
{
  // write number of cell type entries
  vtk_file << "\nCELL_TYPES " << surrogate_.n_pins_ * n_sections_ << "\n";
  // pin loop
  for (size_t pin = 0; pin < surrogate_.n_pins_; pin++) {
    // write the template for each pin
    for (auto v : types_) {
      vtk_file << v << "\n";
    }
  }
  vtk_file << "\n";
} // write_element_types

void SurrogateVtkWriter::write_data(ofstream& vtk_file)
{
  // fuel mesh elements are written first, followed by cladding elements
  // the data needs to be written in a similar matter

  // for each radial section, the data point for that radial ring
  // is repeated azimuthal_res times
  vtk_file << "CELL_DATA " << surrogate_.n_pins_ * n_sections_ << "\n";

  // write temperatures
  if (output_includes_temp_) {
    vtk_file << "SCALARS TEMPERATURE double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";

    // write data for each pin
    for (size_t pin = 0; pin < surrogate_.n_pins_; pin++) {
      if (output_includes_solid_) {
        // write all fuel data first
        for (size_t i = 0; i < n_axial_sections_; i++) {
          for (size_t j = 0; j < n_radial_fuel_sections_; j++) {
            for (size_t k = 0; k < azimuthal_res_; k++) {
              vtk_file << surrogate_.solid_temperature(pin, i, j) << "\n";
            }
          }
        }

        // then write cladding data
        for (size_t i = 0; i < n_axial_sections_; i++) {
          for (size_t j = 0; j < n_radial_clad_sections_; j++) {
            for (size_t k = 0; k < azimuthal_res_; k++) {
              vtk_file << surrogate_.solid_temperature(
                            pin, i, j + n_radial_fuel_sections_)
                       << "\n";
            }
          }
        }
      } // end of solid writing

      // then write fluid data
      if (output_includes_fluid_) {
        for (size_t i = 0; i < n_axial_sections_; ++i) {
          for (size_t j = 0; j < n_fluid_sections_; ++j) {
            vtk_file << surrogate_.fluid_temperature(pin, i) << "\n";
          }
        }
      }
    } // end pin loop
  }

  if (output_includes_density_) {
    vtk_file << "SCALARS DENSITY double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";

    // write data for each pin
    for (size_t pin = 0; pin < surrogate_.n_pins_; pin++) {
      if (output_includes_solid_) {
        // write all fuel data first
        for (size_t i = 0; i < n_axial_sections_; i++) {
          for (size_t j = 0; j < n_radial_fuel_sections_; j++) {
            for (size_t k = 0; k < azimuthal_res_; k++) {
              vtk_file << 0.0 << "\n";
            }
          }
        }

        // then write cladding data
        for (size_t i = 0; i < n_axial_sections_; i++) {
          for (size_t j = 0; j < n_radial_clad_sections_; j++) {
            for (size_t k = 0; k < azimuthal_res_; k++) {
              vtk_file << 0.0 << "\n";
            }
          }
        }
      } // end of solid writing

      // then write fluid data
      if (output_includes_fluid_) {
        for (size_t i = 0; i < n_axial_sections_; ++i) {
          for (size_t j = 0; j < n_fluid_sections_; ++j) {
            vtk_file << surrogate_.fluid_density(pin, i) << "\n";
          }
        }
      }
    } // end pin loop
  }

  if (output_includes_source_) {
    vtk_file << "SCALARS SOURCE double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";

    // Average over the azimuthal sectors for each radial ring.
    // This is consistent with how the source term is used by the surrogate solver.
    xt::xtensor<double, 3> q = xt::mean(surrogate_.source_, 3);

    for (size_t pin = 0; pin < surrogate_.n_pins_; pin++) {
      // write all fuel data first
      if (output_includes_solid_) {
        for (size_t i = 0; i < n_axial_sections_; i++) {
          for (size_t j = 0; j < n_radial_fuel_sections_; j++) {
            for (size_t k = 0; k < azimuthal_res_; k++) {
              vtk_file << q(pin, i, j) << "\n";
            }
          }
        }

        // then write the cladding data
        for (size_t i = 0; i < n_axial_sections_; i++) {
          for (size_t j = 0; j < n_radial_clad_sections_; j++) {
            for (size_t k = 0; k < azimuthal_res_; k++) {
              vtk_file << q(pin, i, j + n_radial_fuel_sections_) << "\n";
            }
          }
        }
      }

      // then write fluid data
      if (output_includes_fluid_) {
        for (size_t i = 0; i < n_axial_sections_; ++i) {
          for (size_t j = 0; j < n_fluid_sections_; ++j) {
            vtk_file << 0.0 << "\n";
          }
        }
      }
    } // end pin loop
  }

} // write_data

xtensor<double, 1> SurrogateVtkWriter::points_for_pin(double x, double y)
{
  // start with the origin-centered template
  xtensor<double, 1> points_out = points_;

  // translate points to pin center
  xt::view(points_out, xt::range(0, _, 3)) += x;
  xt::view(points_out, xt::range(1, _, 3)) += y;

  return points_out;
}

xtensor<double, 1> SurrogateVtkWriter::points()
{
  if (VizRegionType::all == regions_out_) {
    xtensor<double, 3> fuel_pnts = fuel_points();
    xtensor<double, 3> clad_pnts = clad_points();
    xtensor<double, 3> fluid_pnts = fluid_points();
    auto solid_pnts =
      xt::concatenate(xt::xtuple(xt::flatten(fuel_pnts), xt::flatten(clad_pnts)));
    return xt::concatenate(xt::xtuple(solid_pnts, xt::flatten(fluid_pnts)));
  } else if (VizRegionType::solid == regions_out_) {
    xtensor<double, 3> fuel_pnts = fuel_points();
    xtensor<double, 3> clad_pnts = clad_points();
    return xt::concatenate(xt::xtuple(xt::flatten(fuel_pnts), xt::flatten(clad_pnts)));
  } else if (VizRegionType::fluid == regions_out_) {
    return xt::flatten(fluid_points());
  }
}

xtensor<double, 3> SurrogateVtkWriter::fuel_points()
{
  // array to hold all point data
  xt::xarray<double> pnts_out({n_axial_points_, fuel_points_per_plane_, 3}, 0.0);

  // x and y for a single axial plane in the rod
  xt::xarray<double> x = xt::zeros<double>({n_radial_fuel_sections_, azimuthal_res_});
  xt::xarray<double> y = xt::zeros<double>({n_radial_fuel_sections_, azimuthal_res_});

  // generate x/y points for each radial section, beginning with the first non-zero
  // radius, since later we will set the center point tp (0, 0, z) automatically
  for (size_t i = 0; i < n_radial_fuel_sections_; i++) {
    double ring_rad = surrogate_.r_grid_fuel_(i + 1);

    // create a unit ring scaled by radius
    xtensor<double, 2> ring = create_ring(ring_rad, azimuthal_res_);

    // set x/y values for this ring
    xt::view(x, i, xt::all()) = xt::view(ring, 0, xt::all());
    xt::view(y, i, xt::all()) = xt::view(ring, 1, xt::all());
  }

  // flatten radial point arrays
  x = xt::flatten(x);
  y = xt::flatten(y);

  // set the point values for each axial plane for all but the center point, which is
  // always (0, 0, z)
  for (size_t i = 0; i < n_axial_points_; i++) {
    xt::view(pnts_out, i, xt::range(1, _), 0) = x;
    xt::view(pnts_out, i, xt::range(1, _), 1) = y;
    xt::view(pnts_out, i, xt::all(), 2) = surrogate_.z_(i);
  }

  return pnts_out;
}

xtensor<double, 3> SurrogateVtkWriter::clad_points()
{
  // array to hold all point data
  xt::xarray<double> pnts_out({n_axial_points_, clad_points_per_plane_, 3}, 0.0);

  // x and y for a single axial plane in the rod
  xt::xarray<double> x = xt::zeros<double>({n_radial_clad_sections_ + 1, azimuthal_res_});
  xt::xarray<double> y = xt::zeros<double>({n_radial_clad_sections_ + 1, azimuthal_res_});

  // generate x/y points for each radial ring
  for (size_t i = 0; i < n_radial_clad_sections_ + 1; i++) {
    double ring_rad = surrogate_.r_grid_clad_(i);

    // create a unit ring scaled by radius
    xtensor<double, 2> ring = create_ring(ring_rad, azimuthal_res_);

    // set x/y values for this ring
    xt::view(x, i, xt::all()) = xt::view(ring, 0, xt::all());
    xt::view(y, i, xt::all()) = xt::view(ring, 1, xt::all());
  }

  // flatten x,y point arrays
  x = xt::flatten(x);
  y = xt::flatten(y);

  for (size_t i = 0; i < n_axial_points_; i++) {
    // set all points, no center point for cladding
    xt::view(pnts_out, i, xt::all(), 0) = x;
    xt::view(pnts_out, i, xt::all(), 1) = y;
    xt::view(pnts_out, i, xt::all(), 2) = surrogate_.z_[i];
  }

  return pnts_out;
}

xtensor<double, 3> SurrogateVtkWriter::fluid_points()
{
  // array to hold all point data
  xt::xarray<double> pnts_out({n_axial_points_, fluid_points_per_plane_, 3}, 0.0);

  // x and y for a single axial plane in the fluid
  xt::xarray<double> x = xt::zeros<double>({fluid_points_per_plane_});
  xt::xarray<double> y = xt::zeros<double>({fluid_points_per_plane_});

  // create ring on clad surface to set points on clad surface
  xtensor<double, 2> ring = create_ring(surrogate_.clad_outer_radius(), azimuthal_res_);

  xt::view(x, xt::range(0, azimuthal_res_)) = xt::view(ring, 0, xt::all());
  xt::view(y, xt::range(0, azimuthal_res_)) = xt::view(ring, 1, xt::all());

  // set remaining points on boundary
  double half_pitch = surrogate_.pin_pitch() / 2.0;
  for (size_t i = 0; i < 8; ++i) {
    size_t index = i + azimuthal_res_;

    if (i >= 1 && i <= 3) {
      y(index) = half_pitch;
    } else if (i >= 5 && i <= 7) {
      y(index) = -half_pitch;
    }

    if (i >= 3 && i <= 5) {
      x(index) = -half_pitch;
    } else if (i == 0 || i == 1 || i == 7) {
      x(index) = half_pitch;
    }
  }

  for (size_t i = 0; i < n_axial_points_; i++) {
    xt::view(pnts_out, i, xt::all(), 0) = x;
    xt::view(pnts_out, i, xt::all(), 1) = y;
    xt::view(pnts_out, i, xt::all(), 2) = surrogate_.z_[i];
  }

  return pnts_out;
}

xtensor<int, 1> SurrogateVtkWriter::conn()
{
  if (VizRegionType::all == regions_out_) {
    // get both sets of points
    xtensor<int, 4> f_conn = fuel_conn();
    xtensor<int, 4> c_conn = clad_conn();

    // adjust the cladding connectivity by the number of points in the fuel mesh
    xt::view(c_conn, xt::all(), xt::all(), xt::all(), xt::range(1, _)) +=
      fuel_points_per_plane_ * n_axial_points_;

    xtensor<int, 3> fl_conn = fluid_conn();

    // adjust the fluid connectivity by the number of points in the solid mesh,
    // then we need to again correct for the invalid connectivities
    xt::view(fl_conn, xt::all(), xt::all(), xt::range(1, _)) +=
      (fuel_points_per_plane_ + clad_points_per_plane_) * n_axial_points_;
    xt::view(fl_conn, xt::all(), xt::all(), xt::range(7, _)) = INVALID_CONN_;

    auto solid_conn =
      xt::concatenate(xt::xtuple(xt::flatten(f_conn), xt::flatten(c_conn)));
    return xt::concatenate(xt::xtuple(solid_conn, xt::flatten(fl_conn)));

  } else if (VizRegionType::solid == regions_out_) {
    // get both sets of points
    xtensor<int, 4> f_conn = fuel_conn();
    xtensor<int, 4> c_conn = clad_conn();

    // adjust the cladding connectivity by the number of points in the fuel mesh
    xt::view(c_conn, xt::all(), xt::all(), xt::all(), xt::range(1, _)) +=
      fuel_points_per_plane_ * n_axial_points_;

    return xt::concatenate(xt::xtuple(xt::flatten(f_conn), xt::flatten(c_conn)));
  } else if (VizRegionType::fluid == regions_out_) {
    return xt::flatten(fluid_conn());
  }
}

xtensor<int, 1> SurrogateVtkWriter::conn_for_pin(size_t offset)
{
  xt::xarray<int> conn_out = conn_;
  conn_out.reshape({n_sections_, CONN_STRIDE_});

  // get locations of all values less than 0
  xt::xarray<bool> mask = conn_out < 0;

  xt::view(conn_out, xt::all(), xt::range(1, _)) += offset;

  conn_out = xt::flatten(conn_out);
  // no masked_view in this version of xtensor,
  // making do with this for now
  for (size_t i = 0; i < mask.size(); i++) {
    if (mask(i)) {
      conn_out(i) = INVALID_CONN_;
    }
  }
  return conn_out;
}

xtensor<int, 4> SurrogateVtkWriter::fuel_conn()
{
  // size output array
  xtensor<int, 4> cells_out = xt::zeros<int>(
    {n_axial_sections_, n_radial_fuel_sections_, azimuthal_res_, HEX_SIZE_ + 1});

  // generate a base layer to be extended in Z
  xtensor<int, 3> base =
    xt::zeros<int>({n_radial_fuel_sections_, azimuthal_res_, HEX_SIZE_});

  // innermost ring (wedges only); each element is written as a hex element even though
  // it is a wedge element because all other elements except for the innermost ring are
  // hex
  xtensor<int, 2> inner_base = xt::zeros<int>({azimuthal_res_, HEX_SIZE_});

  xt::view(inner_base, xt::all(), 1) = xt::arange(size_t(1), azimuthal_res_ + 1);
  xt::view(inner_base, xt::all(), 2) = xt::arange(size_t(2), azimuthal_res_ + 2);

  // adjust last point id for perioic condition
  xt::view(inner_base, azimuthal_res_ - 1, 2) = 1;

  // copy connectivity of first z-layer to the second and shift
  // by the number of points in a plane
  xt::view(inner_base, xt::all(), xt::range(3, 6)) =
    xt::view(inner_base, xt::all(), xt::range(0, 3));
  xt::strided_view(inner_base, {xt::all(), xt::range(3, 6)}) += fuel_points_per_plane_;

  // set values for the innermost ring
  xt::view(base, 0, xt::all(), xt::all()) = inner_base;

  // outer rings (hexes); create a radial base starting at point zero with a shift
  // between the first and second layer of connectivity equal to the number of
  // fuel points in a plane
  xtensor<int, 2> radial_base = hex_ring(1, azimuthal_res_, fuel_points_per_plane_);

  // set the other rings by shifting the initial base
  // by azimuthal_res_ for each ring
  for (size_t i = 1; i < n_radial_fuel_sections_; i++) {
    xt::view(base, i, xt::all(), xt::all()) = radial_base;
    radial_base += azimuthal_res_;
  }

  // set all axial divs using the base connectivity
  // and shifting by the number of points in a plane
  for (size_t j = 0; j < n_axial_sections_; j++) {
    // set layer and increment connectivity by number of points in axial div
    xt::view(cells_out, j, xt::all(), xt::all(), xt::range(1, _)) = base;
    base += fuel_points_per_plane_;
  }

  // innermost ring is always wedges
  xt::view(cells_out, xt::all(), 0, xt::all(), 0) = WEDGE_SIZE_;

  // the rest are hexes
  xt::view(cells_out, xt::all(), xt::range(1, _), xt::all(), 0) = HEX_SIZE_;

  // first ring should be wedges only, invalidate last two entries
  xt::view(cells_out, xt::all(), 0, xt::all(), xt::range(7, _)) = INVALID_CONN_;

  return cells_out;
}

xtensor<int, 4> SurrogateVtkWriter::clad_conn()
{
  // size output array
  xtensor<int, 4> cells_out = xt::zeros<int>(
    {n_axial_sections_, n_radial_clad_sections_, azimuthal_res_, HEX_SIZE_ + 1});

  // base layer to be extended in Z
  xtensor<int, 3> base =
    xt::zeros<int>({n_radial_clad_sections_, azimuthal_res_, HEX_SIZE_});

  // element rings (hexes only); create a radial base starting at point zero
  // with a shift between the first and second layer of connectivity equal to
  // the number of cladding points in a plane
  xtensor<int, 2> radial_base = hex_ring(0, azimuthal_res_, clad_points_per_plane_);

  // set the other rings by shifting the initial base by azimuthal_res_ for each ring
  for (size_t i = 0; i < n_radial_clad_sections_; i++) {
    xt::view(base, i, xt::all(), xt::all()) = radial_base;
    radial_base += azimuthal_res_;
  }

  // set all axial divs using base for the first layer and
  // shift by the number of points in a plane
  for (size_t i = 0; i < n_axial_sections_; i++) {
    // set layer and increment connectivity by number of points in axial div
    xt::view(cells_out, i, xt::all(), xt::all(), xt::range(1, _)) = base;
    base += clad_points_per_plane_;
  }

  // all are hexes
  xt::view(cells_out, xt::all(), xt::all(), xt::all(), 0) = HEX_SIZE_;

  return cells_out;
}

xtensor<int, 3> SurrogateVtkWriter::fluid_conn()
{
  // size output array
  xtensor<int, 3> cells_out =
    xt::zeros<int>({n_axial_sections_, n_fluid_sections_, HEX_SIZE_ + 1});

  // base layer to be extended in Z
  xtensor<int, 2> base = xt::zeros<int>({n_fluid_sections_, HEX_SIZE_});

  size_t azimuthal_pts_per_quad = azimuthal_res_ / 4 + 1;
  size_t n_fluid_sections_per_quad = n_fluid_sections_ / 4;

  for (size_t i = 0; i < chans_per_rod_; ++i) {
    int start_row = i * n_fluid_sections_per_quad;
    int start_ix = i * (azimuthal_pts_per_quad - 1);
    int corner_ix = azimuthal_res_ + 1 + 2 * i;

    // wedges from outer clad to corner
    for (size_t j = 1; j < n_fluid_sections_per_quad - 1; ++j) {
      base(start_row + j, 0) = start_ix + (j - 1);
      base(start_row + j, 1) = start_ix + j;
      base(start_row + j, 2) = corner_ix;
    }

    // first wedge
    base(start_row, 0) = start_ix;
    base(start_row, 1) = corner_ix;
    base(start_row, 2) = corner_ix - 1;

    // last wedge
    base(start_row + azimuthal_pts_per_quad, 0) = start_ix + azimuthal_pts_per_quad - 1;
    base(start_row + azimuthal_pts_per_quad, 1) = corner_ix + 1;
    base(start_row + azimuthal_pts_per_quad, 2) = corner_ix;
  }

  // correct for periodicity
  base(n_fluid_sections_ - 2, 1) = 0;
  base(n_fluid_sections_ - 1, 0) = 0;
  base(n_fluid_sections_ - 1, 1) = azimuthal_res_;

  for (size_t i = 0; i < n_fluid_sections_; ++i) {
    base(i, 3) = base(i, 0) + n_fluid_sections_;
    base(i, 4) = base(i, 1) + n_fluid_sections_;
    base(i, 5) = base(i, 2) + n_fluid_sections_;
  }

  // set all axial divisions using the base for the first layer and shift by
  // the number of points in a plane
  for (size_t i = 0; i < n_axial_sections_; ++i) {
    xt::view(cells_out, i, xt::all(), xt::range(1, _)) = base;
    base += n_fluid_sections_;
  }

  xt::view(cells_out, xt::all(), xt::all(), 0) = WEDGE_SIZE_;
  xt::view(cells_out, xt::all(), xt::all(), xt::range(7, _)) = INVALID_CONN_;

  return cells_out;
}

xtensor<int, 1> SurrogateVtkWriter::types()
{
  if (VizRegionType::all == regions_out_) {
    xt::xtensor<int, 1> ftypes = xt::flatten(fuel_types());
    xt::xtensor<int, 1> ctypes = xt::flatten(clad_types());
    xt::xtensor<int, 1> fltypes = xt::flatten(fluid_types());

    auto solid_types = xt::concatenate(xt::xtuple(ftypes, ctypes));
    return xt::concatenate(xt::xtuple(solid_types, fltypes));
  } else if (VizRegionType::solid == regions_out_) {
    xtensor<int, 3> ftypes = fuel_types();
    xtensor<int, 3> ctypes = clad_types();

    return xt::concatenate(xt::xtuple(xt::flatten(ftypes), xt::flatten(ctypes)));
  } else if (VizRegionType::fluid == regions_out_) {
    return xt::flatten(fluid_types());
  }
}

xtensor<int, 3> SurrogateVtkWriter::fuel_types()
{
  // size the output array
  xtensor<int, 3> types_out =
    xt::zeros<int>({n_axial_sections_, n_radial_fuel_sections_, azimuthal_res_});

  // the inner ring is wedges
  xt::view(types_out, xt::all(), 0, xt::all()) = WEDGE_TYPE_;
  // the rest are hexes
  xt::view(types_out, xt::all(), xt::range(1, _), xt::all()) = HEX_TYPE_;

  return types_out;
}

xtensor<int, 3> SurrogateVtkWriter::clad_types()
{
  // size output array
  xtensor<int, 3> clad_types_out =
    xt::zeros<int>({n_axial_sections_, n_radial_clad_sections_, azimuthal_res_});
  // all elements are hexes
  clad_types_out = xt::full_like(clad_types_out, HEX_TYPE_);

  return clad_types_out;
}

xtensor<int, 2> SurrogateVtkWriter::fluid_types()
{
  xtensor<int, 2> fluid_types_out =
    xt::zeros<int>({n_axial_sections_, n_fluid_sections_});
  fluid_types_out = xt::full_like(fluid_types_out, WEDGE_TYPE_);
  return fluid_types_out;
}

} // namespace enrico
