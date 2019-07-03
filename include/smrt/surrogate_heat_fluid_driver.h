#ifndef SurrogateHeatFluidDriver_h
#define SurrogateHeatFluidDriver_h

#include <memory>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include <gsl/gsl>

#include "Assembly_Model.h"
#include "Single_Pin_Subchannel.h"
#include "Single_Pin_Conduction.h"

namespace enrico {

//! Driver for coupled subchannel/heat conduction simulations. This driver
//! decouples the fluid and solid phase in an explicit manner (i.e. no
//! iteration between the convective heat transfer linking the phases within
//! a single solve), which is sufficient given that this project targets
//! steady-state simulations via fixed point iteration.
class SurrogateHeatFluidDriver {
public:
  //@{
  //! Typedefs
  using SP_Assembly = std::shared_ptr<Assembly_Model>;
  using RCP_PL = Teuchos::RCP<Teuchos::ParameterList>;
  //@}

  //! fluid temperature
  virtual std::vector<double> fluid_temperature() const { return d_pin_temps; }

  //! fluid density
  virtual std::vector<double> density() const { return d_pin_densities; }

  //! solid temperature
  virtual std::vector<double> solid_temperature() const { return d_solid_temps; }

  double pressure_bc_; //! System pressure in [MPa]

private:
  // >>> DATA
  SP_Assembly d_assembly;

  // Note that for channel-centered equations, Nx and Ny are one greater
  //  than number of pins in each direction.
  int d_Nx, d_Ny, d_Nz;

  // Cross sectional area of each channel
  std::vector<double> d_areas;

  // Mass flow rate (kg/s) in each channel
  std::vector<double> d_mdots;

  //! subchannel solver for a single channel
  std::unique_ptr<Single_Pin_Subchannel> d_pin_subchannel;

  //! heat conduction solver for a single rod
  std::unique_ptr<Single_Pin_Conduction> d_pin_conduction;

  //! coolant temperature in [K] for each channel, of total length given by the
  //! product of the number of pins by the number of axial cells
  std::vector<double> d_pin_temps;

  //! coolant density in [g/cm^3] for each channel, of total length given by the
  //! product of the number of pins by the number of axial cells
  std::vector<double> d_pin_densities;

  //! solid temperature in [K]
  std::vector<double> d_solid_temps;

public:
  // Constructor
  SurrogateHeatFluidDriver(SP_Assembly assembly,
                       RCP_PL subchannel_params,
                       RCP_PL conduction_params,
                       const std::vector<double>& dz);

  //! heat source
  std::vector<double> d_pin_powers;

private:
  //! Set up the sizes of solution arrays
  void generate_arrays();

  int channel_index(int ix, int iy) const
  {
    Expects(ix < d_Nx);
    Expects(iy < d_Ny);
    return ix + d_Nx * iy;
  }

  // Solve the fluid flow equations
  void solve_fluid(const std::vector<double>& power);

  // Solve the solid equations
  void solve_solid(const std::vector<double>& power,
    const std::vector<double>& channel_temp, std::vector<double>& fuel_temp)
};

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
#endif // SurrogateHeatFluidDriver_h

//---------------------------------------------------------------------------//
// end of SurrogateHeatFluidDriver.h
//---------------------------------------------------------------------------//
