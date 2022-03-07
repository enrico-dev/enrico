//! \file base_drivers.h
//! Base classes for single- and coupled-physics drivers
#ifndef ENRICO_DRIVERS_H
#define ENRICO_DRIVERS_H

#include "enrico/comm.h"
#include "enrico/mpi_types.h"
#include "enrico/timer.h"

#include <mpi.h>

#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace enrico {

//! Base class for driver that controls a physics solve
class Driver {
public:
  //! Initializes the solver with the given MPI communicator.
  //! \param comm An existing MPI communicator used to initialize the solver
  explicit Driver(MPI_Comm comm)
    : comm_(comm)
    , timer_driver_setup(comm_)
    , timer_init_step(comm_)
    , timer_solve_step(comm_)
    , timer_write_step(comm_)
    , timer_finalize_step(comm_)
  {
#ifdef _OPENMP
#pragma omp parallel default(none) shared(num_threads)
#pragma omp single
    num_threads = omp_get_num_threads();
#endif
  }

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
  bool active() const { return comm_.active(); }

  Comm comm_; //!< The MPI communicator used to run the solver

  Timer timer_driver_setup;  //!< For the setup function(s) of the external solver
  Timer timer_init_step;     //!< For the init_step() member function
  Timer timer_solve_step;    //!< For the solve_step() member function
  Timer timer_write_step;    //!< For the write_step() member function
  Timer timer_finalize_step; //!< For the finalize_step() member function

  //! Number of OpenMP threads
  int num_threads;
};

//! Contains the mean and standard deviation of an uncertain double float
struct UncertainDouble {
  // Constructors
  UncertainDouble() {}
  UncertainDouble(double nom, double sd)
    : mean{nom}
    , std_dev{sd}
  {}

  double mean{}; //!< Mean value
  double std_dev{}; //!< Standard deviation
};

} // namespace enrico

#endif // ENRICO_DRIVERS_H
