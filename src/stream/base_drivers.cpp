#include "base_drivers.h"

namespace stream {

bool HeatFluidsDriver::active() const
{
  return comm_.comm != MPI_COMM_NULL;
}

bool TransportDriver::active() const
{
  return comm_.comm != MPI_COMM_NULL;
}

} // namespace stream
