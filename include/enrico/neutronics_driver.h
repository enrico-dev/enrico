//! \file neutronics_driver.h
//! Base class for single-physics neutronics solver
#ifndef NEUTRONICS_DRIVER_H
#define NEUTRONICS_DRIVER_H

#include "enrico/driver.h"
#include "enrico/message_passing.h"

#include <gsl/gsl>
#include <xtensor/xtensor.hpp>

#include <vector>

namespace enrico {

using CellHandle = gsl::index;

//! Base class for driver that controls a neutronics solve
class NeutronicsDriver : public Driver {
public:
  explicit NeutronicsDriver(MPI_Comm comm)
    : Driver(comm)
  {}

  virtual ~NeutronicsDriver() = default;

  //! Get energy deposition in each material normalized to a given power
  //! \param power User-specified power in [W]
  //! \return Heat source in each material as [W/cm3]
  virtual xt::xtensor<double, 1> heat_source(double power) const = 0;

  virtual std::vector<CellHandle> find(const std::vector<Position>& positions) = 0;
  virtual void set_density(CellHandle cell, double rho) const = 0;
  virtual void set_temperature(CellHandle cell, double T) const = 0;
  virtual double get_density(CellHandle cell) const = 0;
  virtual double get_temperature(CellHandle cell) const = 0;
  virtual double get_volume(CellHandle cell) const = 0;
  virtual bool is_fissionable(CellHandle cell) const = 0;
  virtual std::size_t n_cells() const = 0;

  // TODO: Remove argument
  virtual void create_tallies(std::size_t n) = 0;
};

} // namespace enrico

#endif // NEUTRONICS_DRIVER_H
