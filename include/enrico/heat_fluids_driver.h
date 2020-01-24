//! \file heat_fluids_driver.h
//! Base class for single-physics heat-fluids solver
#ifndef HEAT_FLUIDS_DRIVER_H
#define HEAT_FLUIDS_DRIVER_H

#include "enrico/driver.h"
#include "enrico/message_passing.h"
#include "xtensor/xtensor.hpp"

#include <cstddef> // for size_t

namespace enrico {

//! Base class for driver that controls a heat-fluids solve
class HeatFluidsDriver : public Driver {
public:
  explicit HeatFluidsDriver(MPI_Comm comm, double pressure_bc);

  virtual ~HeatFluidsDriver() = default;

  //! Whether the calling rank has access to the fields returned
  //! by the 'temperature', 'density', and 'fluid_mask' methods
  //! \return Whether calling rank has access
  virtual bool has_coupling_data() const = 0;

  //! Get the temperature in each region
  //! \return Temperature in each region as [K]
  virtual xt::xtensor<double, 1> temperature() const = 0;

  //! Get the density in each region. The interpretation of this density,
  //! i.e. whether it refers to fluid elements or both fluid and solid
  //! elements, is to the discretion of the particular driver.
  //! \return Temperature in each region as [g/cm^3]
  virtual xt::xtensor<double, 1> density() const = 0;

  //! Get the number of local mesh elements
  // TODO: make pure virtual
  virtual int n_local_elem() const { return 0; }

  //! Get the number of global mesh elements
  // TODO: make pure virtual
  virtual std::size_t n_global_elem() const { return 0; }

  template<typename T>
  std::vector<T> gather(const std::vector<T>& local_field);

  double pressure_bc_; //! System pressure in [MPa]

  //! The displacements of local elements, relative to rank 0. Used in an MPI
  //! Gatherv operation.
  // TODO: Move to private
  std::vector<int32_t> local_displs_;

  //! The number of local elements in each rank.
  // TODO: Move to private
  std::vector<int32_t> local_counts_;

protected:
  //! Initialize the counts and displacements of local elements for each MPI Rank.
  void init_displs();
};

template<typename T>
std::vector<T> HeatFluidsDriver::gather(const std::vector<T>& local_field)
{
  std::vector<T> global_field;

  if (this->active()) {
    // Gather all the local quantities on to the root process
    global_field.resize(this->n_global_elem());
    comm_.Gatherv(local_field.data(),
                  local_field.size(),
                  get_mpi_type<T>(),
                  global_field.data(),
                  local_counts_.data(),
                  local_displs_.data(),
                  get_mpi_type<T>());
  }

  return global_field;
}

} // namespace enrico

#endif // HEAT_FLUIDS_DRIVER_H
