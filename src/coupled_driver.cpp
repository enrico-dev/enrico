#include "enrico/coupled_driver.h"
#include "enrico/driver.h"
#include <gsl/gsl>

namespace enrico {

CoupledDriver::CoupledDriver(MPI_Comm comm, pugi::xml_node node) :
  comm_(comm)
{
  // get required coupling parameters
  power_ = node.child("power").text().as_double();
  max_timesteps_ = node.child("max_timesteps").text().as_int();
  max_picard_iter_ = node.child("max_picard_iter").text().as_int();

  // get optional coupling parameters, and set defaults if not provided
  if (node.child("epsilon"))
    epsilon_ = node.child("epsilon").text().as_double();

  Expects(power_ > 0);
  Expects(max_timesteps_ >= 0);
  Expects(max_picard_iter_ >= 0);
  Expects(epsilon_ > 0);
}

void CoupledDriver::execute()
{
  auto& neutronics = getNeutronicsDriver();
  auto& heat = getHeatDriver();

  // loop over time steps
  for (int i_timestep = 0; i_timestep < max_timesteps_; ++i_timestep)
  {
    std::string msg = "i_timestep: " + std::to_string(i_timestep);
    comm_.message(msg);

    // loop over picard iterations
    for (int i_picard = 0; i_picard < max_picard_iter_; ++i_picard)
    {
      std::string msg = "i_picard: " + std::to_string(i_picard);
      comm_.message(msg);

      if (neutronics.active())
      {
        neutronics.init_step();
        neutronics.solve_step();
        neutronics.write_step(i_timestep, i_picard);
        neutronics.finalize_step();
      }

      comm_.Barrier();

      // Update heat source
      update_heat_source();

      if (heat.active())
      {
        heat.init_step();
        heat.solve_step();
        heat.write_step(i_timestep, i_picard);
        heat.finalize_step();
      }

      comm_.Barrier();

      // Update temperature and density
      update_temperature();
      update_density();

      if (is_converged())
      {
        std::string msg = "converged at i_picard = " + std::to_string(i_picard);
        comm_.message(msg);
        break;
      }
    }
    comm_.Barrier();
  }
  heat.write_step();
}

} // namespace enrico
