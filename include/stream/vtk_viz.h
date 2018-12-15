#include <string>
#include <iostream>
#include <fstream>

#include "stream/heat_driver.h"

namespace stream {

class SurrogateToVtk {

private:
  SurrogateToVtk(const SurrogateHeatDriver *surrogate_ptr);

  const SurrogateHeatDriver* sgate;
  int radial_res;
  // axial
  int n_axial_sections;
  // radial
  int n_radial_fuel_sections, n_radial_clad_sections, n_radial_sections;
  // fuel
  int n_fuel_points;
  // clad
  int n_clad_points;
  // per plane info
  int n_sections_per_plane;
  int fuel_points_per_plane, clad_points_per_plane, points_per_plane;
  // mesh element
  int n_fuel_elements, n_clad_elements, n_mesh_elements;
  // vtk data entry
  int n_fuel_entries_per_plane, n_clad_entries_per_plane, n_entries_per_plane;
  int n_entries;

  // totals
  int total_points;

  public:
  void write_vtk();

  xt::xtensor<double, 3> fuel_points();
};

} // stream
