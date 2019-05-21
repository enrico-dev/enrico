#ifndef Single_Pin_Conduction_h
#define Single_Pin_Conduction_h

#include <vector>

// Trilinos includes
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialDenseVector.hpp"

// SCALE includes
#include "Nemesis/harness/DBC.hh"

// vendored includes
#include <gsl/gsl>

namespace enrico {

//===========================================================================//
/*!
 * \class Single_Pin_Conduction
 * \brief Solve thermal conduction in single fuel pin.
 *
 * This class implements a simple 1-D radial thermal conduction equation
 * for each axial level of a fuel pin.  Axial conduction is ignored and there
 * is currently no gap conductance model.  This model should only be used for
 * very approximate qualitative behavior.
 */
//===========================================================================//

class Single_Pin_Conduction {
public:
  //@{
  //! Typedefs
  using RCP_PL = Teuchos::RCP<Teuchos::ParameterList>;
  using Matrix = Teuchos::SerialDenseMatrix<int, double>;
  using Vector = Teuchos::SerialDenseVector<int, double>;
  using Solver = Teuchos::SerialDenseSolver<int, double>;
  //@}

private:
  // >>> DATA

  // Geometry
  double d_fuel_radius;
  double d_clad_radius;
  double d_delta_r_fuel;
  double d_delta_r_clad;
  std::vector<double> d_delta_z;

  // Thermal conductivities
  double d_k_fuel;
  double d_k_clad;

public:
  // Constructor
  Single_Pin_Conduction(RCP_PL& parameters, const std::vector<double>& delta_z);

  // Set fuel radius (cm)
  void set_fuel_radius(double r)
  {
    Expects(r > 0.0);
    Expects(r < 2.0);
    d_fuel_radius = r;
  }

  // Set clad radius (cm)
  void set_clad_radius(double r)
  {
    Expects(r > 0.0);
    Expects(r < 2.0);
    d_clad_radius = r;
  }

  // Solve conduction equation at each axial level
  void solve(const std::vector<double>& power,
             const std::vector<double>& channel_temperature,
             std::vector<double>& fuel_temperature);
};

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
#endif // Single_Pin_Conduction_h

//---------------------------------------------------------------------------//
// end of Single_Pin_Conduction.h
//---------------------------------------------------------------------------//
