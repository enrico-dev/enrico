#include "enrico/heat_fluids_driver.h"

#include <gsl/gsl>
#include <pugixml.hpp>
#include <xtensor/xadapt.hpp>

namespace enrico {

HeatFluidsDriver::HeatFluidsDriver(MPI_Comm comm, pugi::xml_node node)
  : Driver(comm)
{
  pressure_bc_ = node.child("pressure_bc").text().as_double();
  Expects(pressure_bc_ > 0.0);
}

}
