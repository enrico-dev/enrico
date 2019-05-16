//! \file heat_fluids_driver.h
//! Base class for single-physics heat-fluids solver
#ifndef HEAT_FLUIDS_DRIVER_H
#define HEAT_FLUIDS_DRIVER_H

#include "driver.h"
#include "xtensor/xtensor.hpp"

namespace enrico {

//! Base class for driver that controls a heat-fluids solve
class HeatFluidsDriver : public Driver {
public:
  explicit HeatFluidsDriver(MPI_Comm comm) : Driver(comm){};

  //! Get the temperature in each element
  //! \return Temperature in each element as [K]
  virtual xt::xtensor<double, 1> temperature() const = 0;

  virtual xt::xtensor<double, 1> density() const = 0;

  virtual ~HeatFluidsDriver() = default;
};

} // namespace enrico

#endif // HEAT_FLUIDS_DRIVER_H
