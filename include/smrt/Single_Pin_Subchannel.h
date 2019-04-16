#ifndef Single_Pin_Subchannel_h
#define Single_Pin_Subchannel_h

#include <vector>

// Trilinos includes
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

// SCALE includes
#include "Nemesis/harness/DBC.hh"

namespace enrico {

//===========================================================================//
/*!
 * \class Single_Pin_Subchannel
 * \brief Solve subchannel flow in single channel.
 *
 * This class implements a two-equation subchannel model involving
 * conservation of energy and axial momentum.  No lateral flow is accounted
 * for.  The model excludes friction momentum losses and any spacer grid
 * effects.  It is useful for giving qualitatively-correct behavior but
 * should not be used for actual analysis.
 */
//===========================================================================//
class Single_Pin_Subchannel {
public:
  //@{
  //! Typedefs
  using RCP_PL = Teuchos::RCP<Teuchos::ParameterList>;
  //@}

  enum Verbosity { NONE, LOW, HIGH };

private:
  // >>> DATA

  // Solve parameters
  double d_tol;
  int d_max_iters;
  Verbosity d_verbosity;

  // Axial grid
  std::vector<double> d_delta_z;

  // Channel geometry
  double d_area;

  // Flow conditions
  double d_mdot;
  double d_T_inlet;
  double d_p_exit;

public:
  // Constructor
  Single_Pin_Subchannel(RCP_PL& parameters, const std::vector<double>& delta_z);

  // Set channel area (cm^2)
  void set_channel_area(double area)
  {
    // Channel area internally is needed in m^2
    Require(area > 0.0);
    d_area = 1e-4 * area;
  }

  // Set mass flow rate (kg/s)
  void set_mass_flow_rate(double mdot)
  {
    Require(mdot > 0.0);
    Require(mdot < 2.0);
    d_mdot = mdot;
  }

  // Set inlet temperature (K)
  void set_inlet_temperature(double T)
  {
    Require(T > 0);
    Require(T < 1000);
    d_T_inlet = T;
  }

  // Set exit pressure (Pa)
  void set_exit_pressure(double p)
  {
    Require(p > 0);
    Require(p < 2.2e7);
    d_p_exit = p;
  }

  // Solve for single subchannel given power distribution
  void solve(const std::vector<double>& power,
             std::vector<double>& temperature,
             std::vector<double>& density);
};

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
#endif // Single_Pin_Subchannel_h

//---------------------------------------------------------------------------//
// end of Single_Pin_Subchannel.h
//---------------------------------------------------------------------------//
