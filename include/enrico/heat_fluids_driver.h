//! \file heat_fluids_driver.h
//! Base class for single-physics heat-fluids solver
#ifndef HEAT_FLUIDS_DRIVER_H
#define HEAT_FLUIDS_DRIVER_H

#include "enrico/driver.h"
#include "xtensor/xtensor.hpp"

namespace enrico {

//! Base class for driver that controls a heat-fluids solve
class HeatFluidsDriver : public Driver {
public:
  explicit HeatFluidsDriver(MPI_Comm comm)
    : Driver(comm){};

  virtual ~HeatFluidsDriver() = default;

  //! Get the temperature in each region
  //! \return Temperature in each region as [K]
  virtual xt::xtensor<double, 1> temperature() const = 0;

  //! Get the density in each region
  //! \return Temperature in each region as [g/cm^3]
  virtual xt::xtensor<double, 1> density() const = 0;

  //! States whether each region is in fluid
  //! \return For each region, 1 if region is in fluid and 0 otherwise
  virtual xt::xtensor<int, 1> fluid_mask() const = 0;
};

} // namespace enrico

#endif // HEAT_FLUIDS_DRIVER_H
