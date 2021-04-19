#include "enrico/coupled_driver.h"

#include "enrico/comm_split.h"
#include "enrico/driver.h"
#include "enrico/error.h"
#include "iapws/iapws.h"

#ifdef USE_NEK5000
#include "enrico/nek5000_driver.h"
#endif

#ifdef USE_NEKRS
#include "enrico/nekrs_driver.h"
#endif

#include "enrico/openmc_driver.h"
#ifdef USE_SHIFT
#include "enrico/shift_driver.h"
#endif
#include "enrico/surrogate_heat_driver.h"

#include <gsl/gsl>
#include <xtensor/xbuilder.hpp> // for empty
#include <xtensor/xnorm.hpp>    // for norm_l1, norm_l2, norm_linf

#include <algorithm> // for copy
#include <iomanip>
#include <memory> // for make_unique
#include <string>

// For gethostname
#ifdef _WIN32
#include <winsock.h>
#else
#include <unistd.h>
#endif

namespace enrico {

CoupledDriver::CoupledDriver(MPI_Comm comm, pugi::xml_node node)
  : comm_(comm)
{
  auto neut_node = node.child("neutronics");
  auto heat_node = node.child("heat_fluids");
  auto coup_node = node.child("coupling");

  // get required coupling parameters
  power_ = coup_node.child("power").text().as_double();
  max_timesteps_ = coup_node.child("max_timesteps").text().as_int();
  max_picard_iter_ = coup_node.child("max_picard_iter").text().as_int();

  // get optional coupling parameters, using defaults if not provided
  if (coup_node.child("epsilon"))
    epsilon_ = coup_node.child("epsilon").text().as_double();

  // Determine relaxation parameters for heat source, temperature, and density
  auto set_alpha = [](pugi::xml_node node, double& alpha) {
    if (node) {
      std::string s = node.child_value();
      if (s == "robbins-monro") {
        alpha = ROBBINS_MONRO;
      } else {
        alpha = node.text().as_double();
        Expects(alpha > 0 && alpha <= 1.0);
      }
    }
  };
  set_alpha(coup_node.child("alpha"), alpha_);
  set_alpha(coup_node.child("alpha_T"), alpha_T_);
  set_alpha(coup_node.child("alpha_rho"), alpha_rho_);

  // check for convergence norm
  if (coup_node.child("convergence_norm")) {
    std::string s = coup_node.child_value("convergence_norm");
    if (s == "L1") {
      norm_ = Norm::L1;
    } else if (s == "L2") {
      norm_ = Norm::L2;
    } else if (s == "Linf") {
      norm_ = Norm::LINF;
    } else {
      throw std::runtime_error{"Invalid value for <convergence_norm>"};
    }
  }

  if (coup_node.child("temperature_ic")) {
    std::string s = coup_node.child_value("temperature_ic");

    if (s == "neutronics") {
      temperature_ic_ = Initial::neutronics;
    } else if (s == "heat_fluids") {
      temperature_ic_ = Initial::heat;
    } else {
      throw std::runtime_error{"Invalid value for <temperature_ic>"};
    }
  }

  Expects(power_ > 0);
  Expects(max_timesteps_ >= 0);
  Expects(max_picard_iter_ >= 0);
  Expects(epsilon_ > 0);

  // Create communicators
  std::array<int, 2> nodes{neut_node.child("nodes").text().as_int(),
                           heat_node.child("nodes").text().as_int()};
  std::array<int, 2> procs_per_node{neut_node.child("procs_per_node").text().as_int(),
                                    heat_node.child("procs_per_node").text().as_int()};
  std::array<Comm, 2> driver_comms;
  Comm intranode_comm; // Not used in current comm scheme
  Comm coupling_comm;  // Not used in current comm scheme

  get_driver_comms(
    comm_, nodes, procs_per_node, driver_comms, intranode_comm, coupling_comm);

  auto neutronics_comm = driver_comms[0];
  auto heat_comm = driver_comms[1];

  // Instantiate neutronics driver
  std::string neut_driver = neut_node.child_value("driver");
  if (neut_driver == "openmc") {
    neutronics_driver_ = std::make_unique<OpenmcDriver>(neutronics_comm.comm);
  } else if (neut_driver == "shift") {
#ifdef USE_SHIFT
    neutronics_driver_ = std::make_unique<ShiftDriver>(comm, neut_node);
#else
    throw std::runtime_error{"ENRICO has not been built with Shift support enabled."};
#endif
  } else {
    throw std::runtime_error{"Invalid value for <neutronics><driver>"};
  }

  // Instantiate heat-fluids driver
  std::string s = heat_node.child_value("driver");
  if (s == "nek5000") {
#ifdef USE_NEK5000
    heat_fluids_driver_ = std::make_unique<Nek5000Driver>(heat_comm.comm, heat_node);
#else
    throw std::runtime_error{
      "nek5000 was specified as a solver, but is not enabled in this build of ENRICO"};
#endif
  } else if (s == "nekrs") {
#ifdef USE_NEKRS
    heat_fluids_driver_ = std::make_unique<NekRSDriver>(heat_comm.comm, heat_node);
#else
    throw std::runtime_error{
      "nekrs was specified as a solver, but is not enabled in this build of ENRICO"};
#endif
  } else if (s == "surrogate") {
    heat_fluids_driver_ =
      std::make_unique<SurrogateHeatDriver>(heat_comm.comm, heat_node);
  } else {
    throw std::runtime_error{"Invalid value for <heat_fluids><driver>"};
  }

  // Send rank of neutronics root to all procs
  neutronics_root_ = this->get_neutronics_driver().comm_.is_root() ? comm_.rank : -1;
  MPI_Allreduce(MPI_IN_PLACE, &neutronics_root_, 1, MPI_INT, MPI_MAX, comm_.comm);

  // Send rank of heat root to all procs
  heat_root_ = this->get_heat_driver().comm_.is_root() ? comm_.rank : -1;
  MPI_Allreduce(MPI_IN_PLACE, &heat_root_, 1, MPI_INT, MPI_MAX, comm_.comm);

  // Send number of global elements to all procs
  n_global_elem_ = this->get_heat_driver().n_global_elem();
  comm_.broadcast(n_global_elem_, heat_root_);

  comm_report();

  init_mappings();
  init_tallies();
  init_volumes();

  // elem_fluid_mask_ must be initialized before cell_fluid_mask_!
  init_elem_fluid_mask();
  init_cell_fluid_mask();

  init_temperatures();
  init_heat_source();
}

void CoupledDriver::execute()
{
  auto& neutronics = get_neutronics_driver();
  auto& heat = get_heat_driver();

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

      // Update heat source.
      // On the first iteration, there is no previous iterate of heat source,
      // so we can't apply underrelaxation at that point
      update_heat_source(i_timestep_ > 0 || i_picard_ > 0);

      if (heat.active()) {
        heat.init_step();
        heat.solve_step();
        heat.write_step(i_timestep_, i_picard_);
        heat.finalize_step();
      }

      comm_.Barrier();

      // Update temperature and density
      // At this point, there is always a previous iterate of temperature and density
      // (as assured by the initial conditions set in init_temperature and init_density)
      // so we always apply underrelaxation here.
      update_temperature_and_density(true);

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

double CoupledDriver::temperature_norm(Norm norm)
{
  switch (norm) {
  case Norm::L1: {
    return xt::norm_l1(temperatures_ - temperatures_prev_)();
  }
  case Norm::L2: {
    return xt::norm_l2(temperatures_ - temperatures_prev_)();
  }
  default: {
    return xt::norm_linf(temperatures_ - temperatures_prev_)();
  }
  }
}

bool CoupledDriver::is_converged()
{
  bool converged;
  double norm;

  // The heat root has the global temperature data
  if (comm_.rank == heat_root_) {
    norm = this->temperature_norm(norm_);
    converged = norm < epsilon_;
  }

  comm_.broadcast(converged, heat_root_);
  comm_.broadcast(norm, heat_root_);

  std::string msg = "temperature norm: " + std::to_string(norm);
  comm_.message(msg);
  return converged;
}

void CoupledDriver::update_heat_source(bool relax)
{
  auto& neutronics = this->get_neutronics_driver();
  auto& heat = this->get_heat_driver();

  // ****************************************************************************
  // Gather heat source on neutronics root and apply underrelaxation
  // ****************************************************************************

  if (relax && comm_.rank == neutronics_root_) {
    // Store previous heat source solution if more than one iteration has been performed
    // (otherwise there is not an initial condition for the heat source)
    std::copy(heat_source_.begin(), heat_source_.end(), heat_source_prev_.begin());
  }

  if (neutronics.active()) {
    heat_source_ = neutronics.heat_source(power_);
  }

  // Compute the next iterate of the heat source
  if (relax && comm_.rank == neutronics_root_) {
    if (alpha_ == ROBBINS_MONRO) {
      int n = i_picard_ + 1;
      heat_source_ = heat_source_ / n + (1. - 1. / n) * heat_source_prev_;
    } else {
      heat_source_ = alpha_ * heat_source_ + (1.0 - alpha_) * heat_source_prev_;
    }
  }

  // ****************************************************************************
  // Send underrelaxed heat source to all heat ranks and set element temperatures
  // ****************************************************************************

  this->comm_.send_and_recv(heat_source_, heat_root_, neutronics_root_);
  heat.comm_.broadcast(heat_source_);

  if (heat.active()) {
    // Determine displacement for this rank
    auto displacement = heat.local_displs_.at(heat.comm_.rank);
    int n_local_elem = heat.n_local_elem();
    // Set heat source in every element
    for (int32_t local_elem = 0; local_elem < n_local_elem; ++local_elem) {
      int32_t global_elem = local_elem + displacement;
      // Get heat source for this element
      CellHandle cell = elem_to_cell_.at(global_elem);
      err_chk(heat.set_heat_source_at(local_elem, heat_source_.at(cell)),
              "Error setting heat source for local element " +
                std::to_string(local_elem));
    }
  }
}

void CoupledDriver::update_temperature_and_density(bool relax)
{
  comm_.message("Updating temperature and density");

  auto& neutronics = this->get_neutronics_driver();
  auto& heat = this->get_heat_driver();

  // *************************************************************************
  // Gather temperature on heat root and apply underrelaxation
  // *************************************************************************

  if (relax && comm_.rank == heat_root_) {
    // Store previous temperature solution; a previous solution will always be present
    // because a temperature IC is set and the neutronics solver runs first
    std::copy(temperatures_.begin(), temperatures_.end(), temperatures_prev_.begin());
  }

  temperatures_ = heat.temperature();

  if (relax && comm_.rank == heat_root_) {
    if (alpha_T_ == ROBBINS_MONRO) {
      int n = i_picard_ + 1;
      temperatures_ = temperatures_ / n + (1. - 1. / n) * temperatures_prev_;
    } else {
      temperatures_ = alpha_T_ * temperatures_ + (1.0 - alpha_T_) * temperatures_prev_;
    }
  }

  // ****************************************************************************
  // Send underrelaxed temperature to all neutron ranks and set cell temperatures
  // ****************************************************************************

  // Broadcast global_element_temperatures onto all the neutronics procs
  comm_.send_and_recv(temperatures_, neutronics_root_, heat_root_);
  neutronics.comm_.broadcast(temperatures_);

  if (neutronics.active()) {
    // For each neutronics cell, volume average temperatures and set
    for (const auto& kv : cell_to_elems_) {
      // Get volume-average temperature for this cell instance
      double average_temp = 0.0;
      double total_vol = 0.0;
      for (int elem : kv.second) {
        double T = temperatures_[elem];
        double V = elem_volumes_[elem];
        average_temp += T * V;
        total_vol += V;
      }
      average_temp /= total_vol;
      Ensures(average_temp > 0.0);

      // Set temperature for cell instance
      CellHandle cell = kv.first;
      neutronics.set_temperature(cell, average_temp);

      // If cell is in fluid, set density for cell instance
      if (cell_fluid_mask_[cell] == 1) {
        double density = 1e-3 / iapws::nu1(heat.pressure_bc_, average_temp);
        neutronics.set_density(cell, density);
      }
    }
  }
}

void CoupledDriver::init_mappings()
{
  comm_.message("Initializing mappings");

  const auto& heat = this->get_heat_driver();
  auto& neutronics = this->get_neutronics_driver();

  // Get centroids from heat driver and send to all neutronics procs
  auto elem_centroids = heat.centroids();
  comm_.send_and_recv(elem_centroids, neutronics_root_, heat_root_);
  neutronics.comm_.broadcast(elem_centroids);

  if (neutronics.active()) {
    // Get cell handle corresponding to each element centroid
    elem_to_cell_ = neutronics.find(elem_centroids);

    // Create a vector of elements for each neutronics cell
    for (int32_t elem = 0; elem < elem_to_cell_.size(); ++elem) {
      auto cell = elem_to_cell_[elem];
      cell_to_elems_[cell].push_back(elem);
    }

    // Determine number of neutronic cell instances
    n_cells_ = cell_to_elems_.size();
  }

  // Send element -> cell instance mapping to all heat procs
  comm_.send_and_recv(elem_to_cell_, heat_root_, neutronics_root_);
  heat.comm_.broadcast(elem_to_cell_);

  // Send number of cell instances to all heat procs
  comm_.send_and_recv(n_cells_, heat_root_, neutronics_root_);
  heat.comm_.broadcast(n_cells_);
}

void CoupledDriver::init_tallies()
{
  comm_.message("Initializing tallies");

  auto& neutronics = this->get_neutronics_driver();
  if (neutronics.active()) {
    neutronics.create_tallies();
  }
}

void CoupledDriver::init_temperatures()
{
  comm_.message("Initializing temperatures");

  const auto& neutronics = this->get_neutronics_driver();
  const auto& heat = this->get_heat_driver();

  if (comm_.rank == heat_root_) {
    temperatures_.resize({static_cast<unsigned long>(n_global_elem_)});
    temperatures_prev_.resize({static_cast<unsigned long>(n_global_elem_)});
  } else if (comm_.rank == neutronics_root_) {
    temperatures_.resize({static_cast<unsigned long>(n_global_elem_)});
  }

  if (temperature_ic_ == Initial::neutronics) {
    if (comm_.rank == neutronics_root_) {
      // Loop over heat-fluids elements and determine temperature based on
      // corresponding neutronics cell. This mapping assumes that each
      // heat-fluids element is fully contained within a neutronic cell, i.e.,
      // heat-fluids elements are not split between multiple neutronics cells.
      for (gsl::index elem = 0; elem < elem_to_cell_.size(); ++elem) {
        auto cell = elem_to_cell_[elem];
        double T = neutronics.get_temperature(cell);
        temperatures_[elem] = T;
      }
    }
    comm_.send_and_recv(temperatures_, heat_root_, neutronics_root_);
  } else if (temperature_ic_ == Initial::heat) {
    // * This sets temperatures_ on the the coupling_root, based on the
    //   temperatures received from the heat solver.
    // * We do not want to apply underrelaxation here (and at this point,
    //   there is no previous iterate of temperature, anyway).
    update_temperature_and_density(false);
  }

  // In both cases, only temperatures_ was set, so we explicitly set temperatures_prev_
  if (comm_.rank == heat_root_) {
    std::copy(temperatures_.begin(), temperatures_.end(), temperatures_prev_.begin());
  }
}

void CoupledDriver::init_volumes()
{
  comm_.message("Initializing volumes");

  const auto& heat = this->get_heat_driver();
  const auto& neutronics = this->get_neutronics_driver();

  // Gather all the element volumes on heat root and send to all neutronics procs
  elem_volumes_ = heat.volumes();
  this->comm_.send_and_recv(elem_volumes_, neutronics_root_, heat_root_);
  neutronics.comm_.broadcast(elem_volumes_);

  // Volume check
  if (neutronics.comm_.is_root()) {
    for (const auto& kv : cell_to_elems_) {
      CellHandle cell = kv.first;
      double v_neutronics = neutronics.get_volume(cell);
      double v_heatfluids = 0.0;
      for (const auto& elem : kv.second) {
        v_heatfluids += elem_volumes_.at(elem);
      }

      // Display volumes
      std::stringstream msg;
      msg << "Cell " << neutronics.cell_label(cell) << ", V = " << v_neutronics
          << " (Neutronics), " << v_heatfluids << " (Heat/Fluids)";
      comm_.message(msg.str());
    }
  }
}

void CoupledDriver::init_elem_fluid_mask()
{
  comm_.message("Initializing element fluid mask");

  const auto& heat = this->get_heat_driver();
  const auto& neutronics = this->get_neutronics_driver();

  // Get fluid mask and send to all neutronics procs
  elem_fluid_mask_ = heat.fluid_mask();
  comm_.send_and_recv(elem_fluid_mask_, neutronics_root_, heat_root_);
  neutronics.comm_.broadcast(elem_fluid_mask_);
}

void CoupledDriver::init_cell_fluid_mask()
{
  comm_.message("Initializing cell fluid mask");

  // Because init_elem_fluid_mask is *expected* to be called first, it's assumed that
  // all neutron procs will have the correct values for elem_fluid_mask_
  if (this->get_neutronics_driver().active()) {
    auto n = cell_to_elems_.size();
    cell_fluid_mask_.resize({n});

    for (const auto& kv : cell_to_elems_) {
      CellHandle cell = kv.first;
      for (const auto& elem : kv.second) {
        if (elem_fluid_mask_[elem] == 1) {
          cell_fluid_mask_[cell] = 1;
          break;
        }
        cell_fluid_mask_[cell] = 0;
      }
    }
  }
}

void CoupledDriver::init_heat_source()
{
  comm_.message("Initializing heat source");

  if (comm_.rank == neutronics_root_ || this->get_heat_driver().active()) {
    heat_source_ = xt::empty<double>({n_cells_});
    heat_source_prev_ = xt::empty<double>({n_cells_});
  }
}

void CoupledDriver::comm_report()
{
  char c[_POSIX_HOST_NAME_MAX];
  gethostname(c, _POSIX_HOST_NAME_MAX);
  std::string hostname{c};

  Comm world(MPI_COMM_WORLD);

  // Padding for fields
  int hostw = std::max(8UL, hostname.size()) + 2;
  int rankw = 7;

  for (int i = 0; i < world.size; ++i) {
    if (world.rank == i) {
      if (i == 0) {
        std::cout << "[ENRICO]: Communicator layout: " << std::endl;
        std::cout << std::left << std::setw(hostw) << "Hostname" << std::right
                  << std::setw(rankw) << "World" << std::right << std::setw(rankw)
                  << "Coup" << std::right << std::setw(rankw) << "Neut" << std::right
                  << std::setw(rankw) << "Heat" << std::endl;
      }
      std::cout << std::left << std::setw(hostw) << hostname << std::right
                << std::setw(rankw) << world.rank << std::right << std::setw(rankw)
                << comm_.rank << std::right << std::setw(rankw)
                << this->get_neutronics_driver().comm_.rank << std::right
                << std::setw(rankw) << this->get_heat_driver().comm_.rank << std::endl;
    }
    MPI_Barrier(world.comm);
  }
}

} // namespace enrico
