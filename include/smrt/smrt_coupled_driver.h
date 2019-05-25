#ifndef SMRT_COUPLED_DRIVER_H
#define SMRT_COUPLED_DRIVER_H

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

  double power_; //!< Power in [W]

  int max_picard_iter_; //!< Maximum number of Picard iterations
};

//---------------------------------------------------------------------------//
} // end namespace enrico

#endif // SMRT_COUPLED_DRIVER_H

