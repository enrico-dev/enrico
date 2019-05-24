#include "enrico/heat_fluids_driver.h"
#include <gsl/gsl>

namespace enrico {

HeatFluidsDriver::HeatFluidsDriver(MPI_Comm comm, double pressure)
  : Driver(comm)
  , pressure_(pressure)
{
  Expects(pressure_ > 0.0);
}

}
