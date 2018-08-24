#include "drivers.h"
#include "mpi.h"
#include "nek_interface.h"
#include "openmc/capi.h"
#include "stream_geom.h"

#include <algorithm> // for max, fill, copy
#include <iterator> // for back_inserter
#include <unordered_set>

namespace stream {

bool HeatFluidsDriver::active() const {
  return comm_.comm != MPI_COMM_NULL;
}

bool TransportDriver::active() const {
  return comm_.comm != MPI_COMM_NULL;
}

} // namespace stream
