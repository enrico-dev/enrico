#ifndef STREAM_OPENMC_DRIVER_H
#define STREAM_OPENMC_DRIVER_H

#include <map>
#include <unordered_map>
#include <vector>
#include "mpi.h"
#include "openmc/capi.h"
#include "comm.h"
#include "stream_geom.h"
#include "openmc_interface.h"
#include "drivers.h"

namespace stream {

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

#endif //STREAM_OPENMC_DRIVER_H
