#ifndef Assembly_Model_h
#define Assembly_Model_h

#include <vector>

// Trilinos includes
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

// SCALE includes
#include "Nemesis/harness/DBC.hh"

// vendored includes
#include <gsl/gsl>

namespace enrico {

//===========================================================================//
/*!
 * \class Assembly_Model
 * \brief High-level description of fuel assembly for surrogate models
 */
//===========================================================================//

class Assembly_Model {
public:
  enum PIN_TYPE { FUEL, GUIDE };

private:
  // >>> DATA
  int d_Nx;
  int d_Ny;
  std::vector<PIN_TYPE> d_pin_map;
  std::vector<double> d_x_edges;
  std::vector<double> d_y_edges;
  double d_height;

  // Cylinder radii
  double d_fuel_radius;
  double d_clad_radius;
  double d_guide_radius;

public:
  // Constructor
  Assembly_Model(const std::vector<PIN_TYPE>& pin_map,
                 const std::vector<double>& x_edges,
                 const std::vector<double>& y_edges,
                 double height);

  // Constructor from Teuchos ParameterList
  Assembly_Model(Teuchos::RCP<Teuchos::ParameterList> params);

  // Set fuel pin radius
  void set_fuel_radius(double fr)
  {
    Expects(fr > 0);
    d_fuel_radius = fr;
  }

  // Set clad outer radius
  void set_clad_radius(double cr)
  {
    Expects(cr > 0);
    d_clad_radius = cr;
  }

  // Set guide tube outer radius
  void set_guide_radius(double gr)
  {
    Expects(gr > 0);
    d_guide_radius = gr;
  }

  // Accessors
  double fuel_radius() const { return d_fuel_radius; }
  double clad_radius() const { return d_clad_radius; }
  double guide_radius() const { return d_guide_radius; }
  const std::vector<double>& x_edges() const { return d_x_edges; }
  const std::vector<double>& y_edges() const { return d_y_edges; }
  int num_pins_x() const { return d_Nx; }
  int num_pins_y() const { return d_Ny; }
  int num_pins() const { return d_Nx * d_Ny; }

  // Convert (i,j) to cardinal pin index
  int pin_id(int i, int j) const
  {
    Expects(i < d_Nx);
    Expects(j < d_Ny);
    return i + d_Nx * j;
  }

  // Convert (i,j,k) to cardinal index
  int index(int i, int j, int k) const
  {
    Expects(i < d_Nx);
    Expects(j < d_Ny);
    return pin_id(i, j) + num_pins() * k;
  }

  // Type of pin at (i,j) location
  PIN_TYPE pin_type(int i, int j) const
  {
    Expects(i < d_Nx);
    Expects(j < d_Ny);
    return d_pin_map[pin_id(i, j)];
  }

  // Assembly height
  double height() const { return d_height; }

  // Flow area (cm^2)
  double flow_area(int i, int j) const;
};

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
#endif // Assembly_Model_h

//---------------------------------------------------------------------------//
// end of Assembly_Model.h
//---------------------------------------------------------------------------//
