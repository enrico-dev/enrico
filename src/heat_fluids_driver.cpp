#include "enrico/heat_fluids_driver.h"
#include <gsl/gsl>

namespace enrico {

HeatFluidsDriver::HeatFluidsDriver(MPI_Comm comm, double pressure_bc)
  : Driver(comm)
  , pressure_bc_(pressure_bc)
{
  Expects(pressure_bc_ > 0.0);
}

}
