//! \file neutronics_driver.h
//! Base class for single-physics neutronics solver
#ifndef NEUTRONICS_DRIVER_H
#define NEUTRONICS_DRIVER_H

#include "driver.h"
#include "xtensor/xtensor.hpp"

namespace enrico {

//! Base class for driver that controls a neutronics solve
class NeutronicsDriver : public Driver {
public:
  explicit NeutronicsDriver(MPI_Comm comm)
    : Driver(comm)
  {}

  //! Get energy deposition in each material normalized to a given power
  //! \param power User-specified power in [W]
  //! \return Heat source in each material as [W/cm3]
  virtual xt::xtensor<double, 1> heat_source(double power) const = 0;
};

} // namespace enrico

#endif // NEUTRONICS_DRIVER_H
