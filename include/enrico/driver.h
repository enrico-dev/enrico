//! \file base_drivers.h
//! Base classes for single- and coupled-physics drivers
#ifndef ENRICO_DRIVERS_H
#define ENRICO_DRIVERS_H

#include "enrico/comm.h"
#include "enrico/message_passing.h"

#include <mpi.h>

#include <vector>

namespace enrico {

//! Base class for driver that controls a physics solve
class Driver {
public:
  //! Initializes the solver with the given MPI communicator.
  //! \param comm An existing MPI communicator used to initialize the solver
  explicit Driver(MPI_Comm comm)
    : comm_(comm)
  {}

  //! Performs the necessary initialization for this solver in one Picard iteration
  virtual void init_step() {}

  //! Performs the solve in one Picard iteration
  virtual void solve_step() {}

  //! Write results for physics solve for given timestep and iteration
  //! \param timestep timestep index
  //! \param iteration iteration index
  virtual void write_step(int timestep, int iteration) {}

  //! Write results for a physics solve at the end of the coupled simulation
  void write_step() { this->write_step(-1, -1); }

  //! Performs the necessary finalization for this solver in one Picard iteration
  virtual void finalize_step() {}

  //! Queries whether the comm for this solver is active
  //! \return True if this comm's solver is not MPI_COMM_NULL
  bool active() const;

  //! Broadcast data across ranks, possibly resizing vector
  //! \param values Values to broadcast (significant at rank 0)
  template<typename T>
  void broadcast(std::vector<T>& values) const;

  Comm comm_; //!< The MPI communicator used to run the solver
};

template<typename T>
void Driver::broadcast(std::vector<T>& values) const
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

#endif // ENRICO_DRIVERS_H
