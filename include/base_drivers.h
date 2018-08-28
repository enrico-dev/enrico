#ifndef STREAM_DRIVERS_H
#define STREAM_DRIVERS_H

#include "mpi.h"
#include "comm.h"

namespace stream {

class HeatFluidsDriver {
public:
  Comm comm_;

  explicit HeatFluidsDriver(MPI_Comm comm) : comm_(comm) {};
  HeatFluidsDriver() {};
  virtual ~HeatFluidsDriver() {};

  virtual void init_step() {};
  virtual void solve_step() {};
  virtual void finalize_step() {};
  bool active() const;
};

class TransportDriver {
public:
  Comm comm_;

  explicit TransportDriver(MPI_Comm comm) : comm_(comm) {};
  TransportDriver() {};
  virtual ~TransportDriver() {};

  virtual void init_step() {};
  virtual void solve_step() {};
  virtual void finalize_step() {};
  bool active() const;
};

class CoupledDriver {
public:
  Comm comm_;

  TransportDriver transport_driver_;
  HeatFluidsDriver heat_fluids_driver_;

  explicit CoupledDriver(MPI_Comm coupled_comm, MPI_Comm neutron_comm, MPI_Comm heat_fluids_comm);
  CoupledDriver() {};
  virtual ~CoupledDriver() {};
};

} // namespace stream

#endif //STREAM_DRIVERS_H
