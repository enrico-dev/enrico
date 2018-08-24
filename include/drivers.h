#ifndef STREAM_DRIVERS_H
#define STREAM_DRIVERS_H

#include <map>
#include <unordered_map>
#include <vector>
#include "mpi.h"
#include "openmc/capi.h"
#include "comm.h"
#include "stream_geom.h"
#include "openmc_interface.h"

namespace stream {

// ============================================================================
// Base Classes
// ============================================================================

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

// ============================================================================
// Implementations
// ============================================================================

class OpenmcDriver : public TransportDriver {
public:
  // Constructors and destructors
  OpenmcDriver(int argc, char *argv[], MPI_Comm comm);
  ~OpenmcDriver();

  // Methods
  void init_step();
  void solve_step();
  void finalize_step();
  Position get_mat_centroid(int32_t mat_id) const;

  // Data
  int32_t index_tally_;   //!< Index in tallies array for fission tally
  int32_t index_filter_;  //!< Index in filters arrays for material filter
  std::vector<CellInstance> cells_;
};

} // namespace stream

#endif //STREAM_DRIVERS_H
