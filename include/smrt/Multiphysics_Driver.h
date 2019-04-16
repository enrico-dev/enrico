#ifndef Multiphysics_Driver_h
#define Multiphysics_Driver_h

#include <memory>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Assembly_Model.h"
#include "Multi_Pin_Conduction.h"
#include "Multi_Pin_Subchannel.h"
#include "Neutronics_Solver.h"

namespace enrico {

//===========================================================================//
/*!
 * \class Multiphysics_Driver
 * \brief Simple driver for multiphysics coupling.
 */
/*!
 * \example Driver/test/tstMultiphysics_Driver.cc
 *
 * Test of Multiphysics_Driver.
 */
//===========================================================================//

class Multiphysics_Driver {
public:
  //@{
  //! Typedefs
  using Assembly = enrico::Assembly_Model;
  using SP_Assembly = std::shared_ptr<Assembly>;
  using SP_Subchannel = std::shared_ptr<enrico::Multi_Pin_Subchannel>;
  using SP_Conduction = std::shared_ptr<enrico::Multi_Pin_Conduction>;
  using SP_Neutronics = std::shared_ptr<enrico::Neutronics_Solver>;
  using RCP_PL = Teuchos::RCP<Teuchos::ParameterList>;
  using Vec_Dbl = std::vector<double>;
  //@}

private:
  // >>> DATA
  SP_Assembly d_assembly;
  SP_Subchannel d_subchannel;
  SP_Conduction d_conduction;
  SP_Neutronics d_neutronics;

  int d_Nx, d_Ny, d_Nz, d_N;

  double d_tol;
  int d_max_iters;
  double d_damping;
  std::string d_verbosity;

  double d_power_magnitude;
  Vec_Dbl d_power;
  Vec_Dbl d_fuel_temperature;
  Vec_Dbl d_coolant_temperature;
  Vec_Dbl d_coolant_density;

public:
  // Constructor
  Multiphysics_Driver(SP_Assembly assembly, RCP_PL params, const Vec_Dbl& dz);

  // Solve problem
  void solve();

  // Accessors
  const Vec_Dbl& power() const { return d_power; }
  const Vec_Dbl& fuel_temperature() const { return d_fuel_temperature; }
  const Vec_Dbl& coolant_temperature() const { return d_coolant_temperature; }
  const Vec_Dbl& coolant_density() const { return d_coolant_density; }
};

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
#endif // Multiphysics_Driver_h

//---------------------------------------------------------------------------//
// end of Multiphysics_Driver.h
//---------------------------------------------------------------------------//
