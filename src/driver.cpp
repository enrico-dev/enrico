#include "enrico/driver.h"

namespace enrico {

bool Driver::active() const
{
  return comm_.active();
}

} // namespace enrico
