#include "enrico/base_drivers.h"

namespace enrico {

bool HeatFluidsDriver::active() const
{
  return comm_.comm != MPI_COMM_NULL;
}

bool TransportDriver::active() const
{
  return comm_.comm != MPI_COMM_NULL;
}

} // namespace enrico
