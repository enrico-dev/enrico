//! \file neutronics_driver.h
//! Base class for single-physics neutronics solver
#ifndef NEUTRONICS_DRIVER_H
#define NEUTRONICS_DRIVER_H

#include "enrico/driver.h"
#include "enrico/message_passing.h"

#include "xtensor/xtensor.hpp"

namespace enrico {

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

  //! Broadcast data across ranks, possibly resizing vector
  //! \param values Values to broadcast (significant at rank 0)
  template<typename T>
  void broadcast(std::vector<T>& values) const;
};

template<typename T>
void NeutronicsDriver::broadcast(std::vector<T>& values) const
{
  if (this->active()) {
    // First broadcast the size of the vector
    int n = values.size();
    comm_.Bcast(&n, 1, MPI_INT);

    // Resize vector (for rank != 0) and broacast data
    if (values.size() != n)
      values.resize(n);
    comm_.Bcast(values.data(), n, get_mpi_type<T>());
  }
}

} // namespace enrico

#endif // NEUTRONICS_DRIVER_H
