#include "enrico/coupled_driver.h"
#include "enrico/driver.h"

namespace enrico {

CoupledDriver::CoupledDriver(MPI_Comm comm, pugi::xml_node node) :
  comm_(comm)
{
  // get coupling parameters
  power_ = node.child("power").text().as_double();
  max_timesteps_ = node.child("max_timesteps").text().as_int();
  max_picard_iter_ = node.child("max_picard_iter").text().as_int();

  if (power_ < 0)
    throw std::runtime_error("Power must be non-negative!");

  if (max_timesteps_ < 0)
    throw std::runtime_error("Maximum number of timesteps must be non-negative!");

  if (max_picard_iter_ < 0)
    throw std::runtime_error("Maximum number of Picard iterations must be non-negative!");
}

void CoupledDriver::execute()
{
  auto& neutronics = getNeutronicsDriver();
  auto& heat = getHeatDriver();

  // loop over time steps
  for (int i_timestep = 0; i_timestep < max_timesteps_; ++i_timestep)
  {
    // loop over picard iterations
    for (int i_picard = 0; i_picard < max_picard_iter_; ++i_picard)
    {
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
    }
  }

  heat.write_step();
}

} // namespace enrico
