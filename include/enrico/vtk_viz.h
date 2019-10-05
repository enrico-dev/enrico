#include <fstream>
#include <iostream>
#include <string>

#include "enrico/surrogate_heat_driver.h"

using std::ofstream;
using xt::xtensor;

namespace enrico {

//! Class providing VTK write capabilities for the surrogate T/H solver.
class SurrogateVtkWriter {

  friend SurrogateHeatDriver;

public:
  //! Type of data to output. Valid options are none, source (the kappa
  //! fission distribution), temp, density, and all of the above.
  enum class VizDataType { none = 0, source = 1, temp = 2, density = 3, all = 4 };

  //! Region of space to output data for. Valid options are none, solid
  //! (the fuel and cladding), fluid, and all of the above.
  enum class VizRegionType { none = 0, solid = 1, fluid = 2, all = 3};

public:
  //! Write the surrogate model to VTK
  void write(std::string filename = "magnolia.vtk");

private:
  //! Initializes the surrogate to VTK writer with a surrogate model.
  //! Can only be called within the SurrogateHeatDriver.
  //!
  //! \param surrogate_ptr    Pointer to the surrogate to write
  //! \param t_res            Radial resolution of the generated VTK mesh
  //! \param regions_to_write Description of spatial regions to write
  //! \param data_to_write    Description of solution data to write
  SurrogateVtkWriter(const SurrogateHeatDriver& surrogate_ptr,
                     size_t t_res,
                     const std::string& regions_to_write,
                     const std::string& data_to_write);

  //! Set the number of sections, or cells, that appear in the various
  //! regions of space we may be plotting. When these are described in an
  //! "axial" or "planar" sense, they refer to the number of 2-D regions of
  //! space that would define the faces of polygonal cells in 3-D space.
  //! In other words, if we refer to 8 axial sections and 3 planar sections, the
  //! 3-D space formed by the extrusion of the plaanr sections along the axial
  //! direction would consist of 24 total sections.
  void set_number_of_sections();

  //! Set the number of points, or vertices, that are required to define
  //! the geometry
  void set_number_of_points();

  //! Set the number of entries in the connectivity matrix, which depend on the
  //! mesh element types (wedges, hexes, etc.) and the number of elements
  void set_number_of_entries();

  //! Write a vtk header for an unstructured grid
  void write_header(ofstream& vtk_file);

  //! Write points to the vtk file
  void write_points(ofstream& vtk_file);

  //! Write the wedge/hex element connectivity to the vtk file
  void write_element_connectivity(ofstream& vtk_file);

  //! Write the wedge/hex element types to the vtk file
  void write_element_types(ofstream& vtk_file);

  //! Write requested data to the vtk file
  void write_data(ofstream& vtk_file);

  //! Generate fuel mesh points
  //! \return fuel points (axial, radial_rings, xyz)
  xtensor<double, 3> fuel_points();

  //! Generate cladding mesh points
  //! \return cladding points (axial, radial_rings, xyz)
  xtensor<double, 3> clad_points();

  //! Generate fluid mesh points; the fluid points are ordered in a counterclockwise
  //! manner beginning with the points on the surface of the cladding followed by the
  //! eight points defining the corners of the four subchannels surrounding a pin.
  //! This is the "specified map" used to describe the second indexing in this class.
  //! \return fluid points (axial, specified map, xyz)
  xtensor<double, 3> fluid_points();

  //! Return 1-D array of points for writing (xyz...) (ordered radially, axially)
  //! \return 1-D array of all points in the model
  xtensor<double, 1> points();

  //! Return 1-D array of points, translated to a pin center
  //! \return 1-D array of points translated to a center x,y
  xtensor<double, 1> points_for_pin(double x, double y);

  //! Return connectivity for a pin with an offset
  xtensor<int, 1> conn_for_pin(size_t offset);

  //! Generate fuel connectivity (axial, radial, res, conn)
  //! \return fuel element connectivity (axial, radial, res, conn)
  xtensor<int, 4> fuel_conn();

  //! Generate cladding connectivity
  //! \return cladding element connectivity (axial, radial, res, conn)
  xtensor<int, 4> clad_conn();

  //! Generate fluid connectivity
  //! \return fluid element connectivity (axial, specified map, conn)
  xtensor<int, 3> fluid_conn();

  //! Return all connectivity values for writing
  //! \return 1-D array of connectivity, strided by conn_stride_ per element (ordered
  //! radially, axially)
  xtensor<int, 1> conn();

  //! Generate fuel vtk mesh element (cell) types
  //! \return fuel cell types (axial, radial, azimuthal)
  xtensor<int, 3> fuel_types();

  //! Generate cladding vtk mesh element (cell) types
  //! \return cladding cell types (axial, radial, azimuthal)
  xtensor<int, 3> clad_types();

  //! Generate fluid vtk mesh element (cell) types
  //! \return fluid cell types (axial, specified map, azimuthal)
  xtensor<int, 2> fluid_types();

  //! Return all element type values for writing
  //! \return 1-D array of types, one for each element (ordered planar, axially)
  xtensor<int, 1> types();

  const SurrogateHeatDriver& surrogate_; //!< reference to surrogate
  size_t azimuthal_res_;                 //!< azimuthal resolution
  VizDataType data_out_;                 //!< output region
  VizRegionType regions_out_;            //!< output data

  //! Whether the output region contains the fluid region
  bool output_includes_fluid_;

  //! Whether the output region contains the solid region
  bool output_includes_solid_;

  //! Whether the output data contains temperature
  bool output_includes_temp_;

  //! Whether the output data contains density
  bool output_includes_density_;

  //! Whether the output data contains the fission source
  bool output_includes_source_;

  //! Stride in connectivity, equal to the hex size plus one if writing the solid
  //! regions and equal to the wedge size plus one is writing the fluid regions.
  size_t conn_stride_;

  //! Number of coolant channels surrounding each rod
  const size_t chans_per_rod_ = 4;

  //!< template of xyz values for the mesh, centerd on the origin
  xtensor<double, 1> points_;

  //!< template of mesh element connectivity for a single pin
  xtensor<int, 1> conn_;

  //!< template of mesh element types for a single pin
  xtensor<int, 1> types_;

  //! number of axial sections for a single rod
  size_t n_axial_sections_;

  //! number of radial fuel sections (before division by azimuthal refinement)
  //! on each plane for a single rod
  size_t n_radial_fuel_sections_;

  //! number of radial cladding sections (before division by azimuthal refinement)
  //! on each plane for a single rod
  size_t n_radial_clad_sections_;

  //! number of fluid sections on each plane for a single rod
  size_t n_fluid_sections_;

  //! total number of sections for a single rod, which represents the multiplication
  //! of the number of axial sections by the number of planar sections
  size_t n_sections_;

  //! number of axial points (planes) for a single rod
  size_t n_axial_points_;

  //! number of fuel points on each plane for a single rod
  size_t fuel_points_per_plane_;

  //! number of cladding points on each plane for a single rod
  size_t clad_points_per_plane_;

  //! number of fluid points on each plane
  size_t fluid_points_per_plane_;

  //! total number of points required to represent each pin
  size_t n_points_;

  //!< number of fuel elements in mesh
  size_t n_fuel_elements_;

  //!< number of cladding elements in the mesh
  size_t n_clad_elements_;

  //! number of fuel connectivity entries on each plane in a single rod;
  //! the innermost ring of fuel sections are wedges, while all others are hexes.
  size_t n_fuel_entries_per_plane_;

  //! number of cladding connectivity entries on each plane in a single rod;
  //! all cladding sections are hexes.
  size_t n_clad_entries_per_plane_;

  //! number of fluid connectivity entries on each plane in a single rod;
  //! all fluid sections are wedges.
  size_t n_fluid_entries_per_plane_;

  //! total number of entries in the connectivity matrix; this number represents
  //! the length of the element-size plus one (to indicate the number of points
  //! in the section) summed over all sections
  size_t n_entries_;
};

} // namespace enrico
