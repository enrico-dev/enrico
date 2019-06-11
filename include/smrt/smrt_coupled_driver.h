#ifndef SMRT_COUPLED_DRIVER_H
#define SMRT_COUPLED_DRIVER_H

#include "Neutronics_Solver.h"
#include "enrico/heat_fluids_driver.h"

namespace enrico {

/**
 * Base class for driver that controls a coupled physics solve involving
 * neutronics and thermal-hydraulics physics. This is intended to be a
 * temporary class to aid in incrementally moving the ShiftNekDriver to be
 * a derived class of CoupledDriver, so most of the methods in this class
 * are created to be similar to CoupledDriver.
 */
class SmrtCoupledDriver {
public:
  SmrtCoupledDriver();

  //! Get reference to neutronics driver
  //! \return reference to driver
  Neutronics_Solver& get_neutronics_driver() const = 0;

  //! Get reference to thermal-fluids driver
  //! \return reference to driver
  HeatFluidsDriver& get_heat_driver() const = 0;

  //! Set the heat source in the thermal-hydraulics solver
  virtual void set_heat_source() {};

  //! Update the temperature for the neutronics solver
  virtual void update_temperature() {}

  //! Update the density for the neutronics solver
  virtual void update_density() {}

  double power_; //!< Power in [W]

  int max_picard_iter_; //!< Maximum number of Picard iterations
};

//---------------------------------------------------------------------------//
} // end namespace enrico

#endif // SMRT_COUPLED_DRIVER_H

