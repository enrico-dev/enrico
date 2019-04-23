#include "enrico/coupled_driver.h"
#include "enrico/driver.h"
#include "xtensor/xnorm.hpp"
#include <gsl/gsl>

namespace enrico {

CoupledDriver::CoupledDriver(MPI_Comm comm, pugi::xml_node node)
  : comm_(comm)
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
  for (i_timestep_ = 0; i_timestep_ < max_timesteps_; ++i_timestep_) {
    std::string msg = "i_timestep: " + std::to_string(i_timestep_);
    comm_.message(msg);

    // loop over picard iterations
    for (i_picard_ = 0; i_picard_ < max_picard_iter_; ++i_picard_) {
      std::string msg = "i_picard: " + std::to_string(i_picard_);
      comm_.message(msg);

      if (neutronics.active()) {
        neutronics.init_step();
        neutronics.solve_step();
        neutronics.write_step(i_timestep_, i_picard_);
        neutronics.finalize_step();
      }

      comm_.Barrier();

      // Update heat source
      update_heat_source();

      if (heat.active()) {
        heat.init_step();
        heat.solve_step();
        heat.write_step(i_timestep_, i_picard_);
        heat.finalize_step();
      }

      comm_.Barrier();

      // Update temperature and density
      update_temperature();
      update_density();

      if (is_converged()) {
        std::string msg = "converged at i_picard = " + std::to_string(i_picard_);
        comm_.message(msg);
        break;
      }
    }
    comm_.Barrier();
  }
  heat.write_step();
}

void CoupledDriver::compute_temperature_norm(const CoupledDriver::Norm& n, double& norm, bool& converged)
{
  switch (n) {
    case Norm::L1: {
      auto n_expr = xt::norm_l1(temperatures_ - temperatures_prev_);
      norm = n_expr();
      break;
    }
    case Norm::L2: {
      auto n_expr = xt::norm_l2(temperatures_ - temperatures_prev_);
      norm = n_expr();
      break;
    }
    default: {
      auto n_expr = xt::norm_linf(temperatures_ - temperatures_prev_);
      norm = n_expr();
      break;
    }
  }

  converged = norm < epsilon_;
}

} // namespace enrico
