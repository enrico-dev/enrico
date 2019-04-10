#include "enrico/driver.h"

namespace enrico {

bool Driver::active() const
{
  return comm_.comm != MPI_COMM_NULL;
}

} // namespace enrico
