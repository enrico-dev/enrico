//! \file base_drivers.h
//! Base classes for single- and coupled-physics drivers
#ifndef ENRICO_DRIVERS_H
#define ENRICO_DRIVERS_H

#include "mpi.h"
#include "comm.h"

namespace enrico {

//! Base class for heat/fluid physics driver
class HeatFluidsDriver {
public:

  Comm comm_; //!< The MPI communicator used to run the solver

  //! Initializes the solver with the given MPI communicator
  //!
  //! This should do one-time initialization for the solver, not the initialization that needs to
  //! happen in each Picard iteration.
  //!
  //! \param comm An existing MPI communicator used to initialize the solver
  explicit HeatFluidsDriver(MPI_Comm comm) : comm_(comm) {};

  //! Default constructor
  HeatFluidsDriver() {};

  //! Virtual destructor
  virtual ~HeatFluidsDriver() {};

  //! Performs the necessary initialization for this solver in one Picard iteration
  virtual void init_step() {};

  //! Performs the solve in one Picard iteration
  virtual void solve_step() {};

  //! Performs the necessary finalization for this solver in one Picard iteration
  virtual void finalize_step() {};

  //! Queries whether the comm for this solver is active
  //! \return True if this comm's solver is not MPI_COMM_NULL
  bool active() const;
};

//! Base class for particle transport physics driver
class TransportDriver {
public:
  Comm comm_; //!< The MPI communicator used to run the solver

  //! Initializes the solver with the given MPI communicator
  //!
  //! This should do one-time initialization for the solver, not the initialization that needs to
  //! happen in each Picard iteration.
  //!
  //! \param comm An existing MPI communicator used to initialize the solver
  explicit TransportDriver(MPI_Comm comm) : comm_(comm) {};

  //! Default constructor
  TransportDriver() {};

  //! Virtual destructor
  virtual ~TransportDriver() {};

  //! Performs the necessary initialization for this solver in one Picard iteration
  virtual void init_step() {};

  //! Performs the solve in one Picard iteration
  virtual void solve_step() {};

  //! Performs the necessary finalization for this solver in one Picard iteration
  virtual void finalize_step() {};

  //! Queries whether the comm for this solver is active
  //! \return True if this comm's solver is not MPI_COMM_NULL
  bool active() const;
};

//! Base class for coupled particle transport and heat/fluids physics
class TransportHeatFluidsDriver {
public:
  Comm comm_; //!< The MPI communicator used to run the solver

  TransportDriver transport_driver_; //!< Driver for particle transport physics
  HeatFluidsDriver heat_fluids_driver_; //!< Driver for heat/fluids physics

  //! Initializes this driver and its member drivers from existing MPI communicators
  //! \param coupled_comm An existing comm for this coupled-physics driver
  //! \param neutron_comm An existing comm for the member transport driver
  //! \param heat_fluids_comm  An existing comm for the member heat/fluids driver
  explicit TransportHeatFluidsDriver(MPI_Comm coupled_comm, MPI_Comm neutron_comm, MPI_Comm heat_fluids_comm);

  //! Default constructor
  TransportHeatFluidsDriver() {};

  virtual void update_heat_source() = 0;

  virtual void update_temperature() = 0;

  //! Virtual desctructor
  virtual ~TransportHeatFluidsDriver(){};
};

} // namespace enrico

#endif //ENRICO_DRIVERS_H
