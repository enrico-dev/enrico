#ifndef STREAM_NEK_DRIVER_H
#define STREAM_NEK_DRIVER_H

#include <map>
#include <unordered_map>
#include <vector>
#include "mpi.h"
#include "openmc/capi.h"
#include "comm.h"
#include "stream_geom.h"
#include "openmc_interface.h"
#include "base_drivers.h"
#include "nek_interface.h"

namespace stream {

class NekDriver : public HeatFluidsDriver {
public:
  explicit NekDriver(MPI_Comm comm);
  ~NekDriver();

  void init_step();
  void solve_step();
  void finalize_step();

  Position get_global_elem_centroid(int global_elem) const;

  int lelg_; //!< upper bound on number of mesh elements
  int lelt_; //!< upper bound on number of mesh elements per rank
  int lx1_; //!< polynomial order of the solution
  int nelgt_; //!< total number of mesh elements
  int nelt_; //!< number of local mesh elements
};

} // namespace stream

#endif //STREAM_NEK_DRIVER_H
