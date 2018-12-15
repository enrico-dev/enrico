#include <string>
#include <iostream>
#include <fstream>

#include "stream/heat_driver.h"

namespace stream {

class SurrogateToVtk {

  friend SurrogateHeatDriver;

  enum class VizDataType {
    none = 0,
    source = 1,
    temp   = 2,
    all    = 3
  };

  enum class VizRegionType {
    none = 0,
    fuel = 1,
    clad = 2,
    all  = 3
  };

private:
  //! Initializes the surrogate to VTK writer with a surrogate model. Can only be called withing the SurrogateHeatDriver.
  //!
  //! \param surrogate_ptr Pointer to the surrogate to write
  //! \param t_rad         Radial resolution of the generated VTK mesh
  SurrogateToVtk(const SurrogateHeatDriver *surrogate_ptr,
                 int t_res,
                 std::string regions_to_write,
                 std::string data_to_write);

public:
  //! Write the surrogate model to VTK
  void write_vtk(std::string filename = "magnolia.vtk");

  //! Generate fuel mesh points
  //! \return fuel points (axial, radial_rings, xyz)
  xt::xtensor<double, 3> fuel_points();
  //! Generate cladding mesh points
  //! \return cladding points (axial, radial_rings, xyz)
  xt::xtensor<double, 3> clad_points();
  //! Return 1-D array of points for writing (xyz...) (ordered radially, axially)
  xt::xtensor<double, 1> points();

  //! Generate fuel connectivity (axial, radial, res, conn)
  //! \return fuel element connectivity (axial, radial, res, conn)
  xt::xtensor<int, 4> fuel_conn();
  //! Generate cladding connectivity
  //! \return cladding element connectivity (axial, radial, res, conn)
  xt::xtensor<int, 4> clad_conn();
  //! Return all connectivity values for writing
  //! \return 1-D array of connectivity, strided by CONN_STRIDE_ per element (ordered radially, axially)
  xt::xtensor<int, 1> conn();

  //! Generate fuel vtk mesh element (cell) types
  //! \return fuel cell types (axia, radial, azimuthal)
  xt::xtensor<int, 3> fuel_types();
  //! Generate cladding vtk mesh element (cell) types
  //! \return cladding cell types (axia, radial, azimuthal)
  xt::xtensor<int, 3> clad_types();
  //! Return all element type values for writing
  //! \return 1-D array of types, one for each element (ordered radially, axially)
  xt::xtensor<int, 1> types();

private:
  // internal variables
  const SurrogateHeatDriver* sgate; //!< pointer to surrogate
  int radial_res; //!< radial resolution;

  VizDataType data_out_;
  VizRegionType regions_out_ ;

  /// SECTIONS \\\
  // axial
  int n_axial_sections; //!< number of axial sections
  int n_axial_points; //!< number of axial points (planes)
  // radial
  int n_radial_fuel_sections; //!< number of radial fuel sections
  int n_radial_clad_sections; //!< number of radial cladding sections
  int n_radial_sections; //!< total number of radial sections
  /// POINTS \\\

  int n_fuel_points; //!< number of points neede to represent the fuel regions
  int n_clad_points; //!< number of points neede to represent the cladding regions
  int n_points; //!< total number of points in the vtk representation

  /// PER PLANE \\\

  int n_sections_per_plane; //!< total number of radial sections in a plane
  int fuel_points_per_plane; //!< number of fuel points in a plane
  int clad_points_per_plane; //!< number of cladding points in a plane
  int points_per_plane; //!< total number of points in a plane

  /// MESH ELEMENTS \\\

  int n_fuel_elements; //!< number of fuel elements in mesh
  int n_clad_elements; //!< number of cladding elements in the mesh
  int n_mesh_elements; //!< total number of mesh elements

  /// CONNECTIVITY ENTRIES \\\

  int n_fuel_entries_per_plane; //!< number of fuel connectivity entries in a plane
  int n_clad_entries_per_plane; //!< number of cladding connectivity entries in a plane
  int n_entries_per_plane; //!< total number of connectivity entries in a plane
  int n_entries; //!< total number of connectivity entries in a plane

};


} // stream
