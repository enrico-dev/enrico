#include "enrico/coupled_driver.h"

#include "enrico/comm_split.h"
#include "enrico/driver.h"
#include "enrico/error.h"

#ifdef USE_NEK
#include "enrico/nek_driver.h"
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
#include <set>
#include <string>

// For gethostname
#ifdef _WIN32
#include <winsock.h>
#else
#include <unistd.h>
#endif

#define PLINE                                                                            \
  std::cout << comm_.rank << " : " << __FILE__ << " : " << __LINE__ << std::endl;

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

  if (coup_node.child("density_ic")) {
    std::string s = coup_node.child_value("density_ic");

    if (s == "neutronics") {
      density_ic_ = Initial::neutronics;
    } else if (s == "heat_fluids") {
      density_ic_ = Initial::heat;
    } else {
      throw std::runtime_error{"Invalid value for <density_ic>"};
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
#ifdef USE_NEK
    heat_fluids_driver_ = std::make_unique<NekDriver>(heat_comm.comm, heat_node);
#else
    throw std::runtime_error{
      "nek5000 was specified as a solver, but is not enabled in this build of ENRICO"};
#endif
  } else if (s == "surrogate") {
    heat_fluids_driver_ =
      std::make_unique<SurrogateHeatDriver>(heat_comm.comm, heat_node);
  } else {
    throw std::runtime_error{"Invalid value for <heat_fluids><driver>"};
  }

  // Discover the ranks that are in each comm
  neutronics_ranks_ = gather_subcomm_ranks(comm_, neutronics_comm);
  heat_ranks_ = gather_subcomm_ranks(comm_, heat_comm);

  // Send rank of neutronics root to all procs
  neutronics_root_ = this->get_neutronics_driver().comm_.is_root() ? comm_.rank : -1;
  MPI_Allreduce(MPI_IN_PLACE, &neutronics_root_, 1, MPI_INT, MPI_MAX, comm_.comm);

  // Send rank of heat root to all procs
  heat_root_ = this->get_heat_driver().comm_.is_root() ? comm_.rank : -1;
  MPI_Allreduce(MPI_IN_PLACE, &heat_root_, 1, MPI_INT, MPI_MAX, comm_.comm);

  comm_report();

  init_mappings();
  init_tallies();
  init_volumes();

  init_cell_fluid_mask();

  // Debug:
  decltype(cell_fluid_mask_) fluid_mask_recv;
  decltype(cells_) cells_recv;
  for (const auto& heat_rank : heat_ranks_) {
    comm_.send_and_recv(cells_recv, neutronics_root_, cells_, heat_rank);
    comm_.send_and_recv(fluid_mask_recv, neutronics_root_, cell_fluid_mask_, heat_rank);
    if (comm_.rank == neutronics_root_) {
      for (gsl::index i = 0; i < cells_recv.size(); ++i) {
        std::cout << "    " << get_neutronics_driver().cell_label(cells_recv.at(i))
                  << " : " << fluid_mask_recv.at(i) << std::endl;
      }
    }
  }

  init_temperatures();
  init_densities();
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
      update_temperature(true);
      update_density(true);

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
  auto& heat = this->get_heat_driver();
  double tnorm = 0;

  // Update global temperature norm
  if (heat.active()) {
    switch (norm) {
    case Norm::L1: {
      double l_sum =
        xt::eval(xt::sum(xt::abs(cell_temperatures_ - cell_temperatures_prev_)))[0];
      MPI_Reduce(&l_sum, &tnorm, 1, MPI_DOUBLE, MPI_SUM, 0, heat.comm_.comm);
      break;
    }
    case Norm::L2: {
      double l_sum =
        xt::eval(xt::sum(xt::square(cell_temperatures_ - cell_temperatures_prev_)))[0];
      double g_sum;
      MPI_Reduce(&l_sum, &g_sum, 1, MPI_DOUBLE, MPI_SUM, 0, heat.comm_.comm);
      tnorm = sqrt(g_sum);
      break;
    }
    case Norm::LINF: {
      double l_max =
        xt::eval(xt::amax(xt::abs(cell_temperatures_ - cell_temperatures_prev_)))[0];
      MPI_Reduce(&l_max, &tnorm, 1, MPI_DOUBLE, MPI_MAX, 0, heat.comm_.comm);
      break;
    }
    }
  }

  return tnorm;
}

bool CoupledDriver::is_converged()
{
  bool converged;
  double norm;

  // The heat root has the global temperature data
  norm = this->temperature_norm(norm_);
  if (comm_.rank == heat_root_) {
    converged = norm < epsilon_;
  }

  comm_.broadcast(converged, heat_root_);
  comm_.broadcast(norm, heat_root_);

  std::stringstream msg;
  msg << "temperature norm: " << norm;
  comm_.message(msg.str());
  return converged;
}

void CoupledDriver::update_heat_source(bool relax)
{
  comm_.message("Updating heat source");
  auto& neutronics = this->get_neutronics_driver();
  auto& heat = this->get_heat_driver();

  // ****************************************************************************
  // Gather heat source on neutronics root and apply underrelaxation
  // ****************************************************************************

  if (relax && heat.active()) {
    std::copy(cell_heat_.cbegin(), cell_heat_.cend(), cell_heat_prev_.begin());
  }

  decltype(cells_) cells_recv;
  decltype(cell_heat_) cell_heat_send;
  xt::xtensor<double, 1> all_cell_heat;

  // The neutronics root sends each heat root the cell-averaged heat sources that it needs
  if (comm_.rank == neutronics_root_) {
    all_cell_heat = neutronics.heat_source(power_);
  }
  for (const auto& heat_rank : heat_ranks_) {
    comm_.send_and_recv(cells_recv, neutronics_root_, cells_, heat_rank);
    cell_heat_send.resize({cells_recv.size()});
    if (comm_.rank == neutronics_root_) {
      for (gsl::index i = 0; i < cells_recv.size(); ++i) {
        auto j = neutronics.cell_index(cells_recv.at(i));
        cell_heat_send.at(i) = all_cell_heat.at(j);
      }
    }
    comm_.send_and_recv(cell_heat_, heat_rank, cell_heat_send, neutronics_root_);
  }

  // On heat rank, update the element heat sources
  if (heat.active()) {
    if (relax) {
      if (alpha_ == ROBBINS_MONRO) {
        int n = i_picard_ + 1;
        cell_heat_ = cell_heat_ / n + (1. - 1. / n) * cell_heat_prev_;
      } else {
        cell_heat_ = alpha_ * cell_heat_ + (1.0 - alpha_) * cell_heat_prev_;
      }
    }
    for (gsl::index i = 0; i < cells_.size(); ++i) {
      for (const auto& e : cell_to_elems_.at(cells_.at(i))) {
        heat.set_heat_source_at(e, cell_heat_.at(i));
      }
    }
  }
}

void CoupledDriver::update_temperature(bool relax)
{
  comm_.message("Updating temperature");
  auto& neutronics = this->get_neutronics_driver();
  auto& heat = this->get_heat_driver();

  // Step 1: On heat ranks, assign the current iterate of T to the previous
  if (relax && heat.active()) {
    std::copy(cell_temperatures_.begin(),
              cell_temperatures_.end(),
              cell_temperatures_prev_.begin());
  }

  // Step 2: On heat ranks, compute cell T
  if (heat.active()) {
    auto elem_temperatures = heat.temperature_local();
    for (gsl::index i = 0; i < cells_.size(); ++i) {
      double avg_temp = 0.0;
      double total_vol = 0.0;
      for (const auto& e : cell_to_elems_.at(cells_.at(i))) {
        double T = elem_temperatures.at(e);
        double V = elem_volumes_.at(e);
        avg_temp += T * V;
        total_vol += V;
      }
      avg_temp /= total_vol;
      Ensures(avg_temp > 0.0);
      cell_temperatures_.at(i) = avg_temp;
    }
    // Apply relaxation to local cell temperatures
    if (relax) {
      if (alpha_T_ == ROBBINS_MONRO) {
        int n = i_picard_ + 1;
        cell_temperatures_ =
          cell_temperatures_ / n + (1. - 1. / n) * cell_temperatures_prev_;
      } else {
        cell_temperatures_ =
          alpha_T_ * cell_temperatures_ + (1.0 - alpha_T_) * cell_temperatures_prev_;
      }
    }
  }

  // Step 3: On neutron ranks, accumulate cell volumes from all heat ranks
  std::map<CellHandle, double> T_dot_V;
  decltype(cells_) cells_recv;
  decltype(cell_volumes_) cell_volumes_recv;
  decltype(cell_temperatures_) cell_temperatures_recv;
  for (const auto& heat_rank : heat_ranks_) {
    comm_.send_and_recv(cells_recv, neutronics_root_, cells_, heat_rank);
    neutronics.comm_.broadcast(cells_recv);

    comm_.send_and_recv(
      cell_temperatures_recv, neutronics_root_, cell_temperatures_, heat_rank);
    neutronics.comm_.broadcast(cell_temperatures_recv);

    comm_.send_and_recv(cell_volumes_recv, neutronics_root_, cell_volumes_, heat_rank);
    neutronics.comm_.broadcast(cell_volumes_recv);

    if (neutronics.active()) {
      for (gsl::index i = 0; i < cells_recv.size(); ++i) {
        auto T = cell_temperatures_recv.at(i);
        auto V = cell_volumes_recv.at(i);
        T_dot_V[cells_recv.at(i)] += T * V;
      }
    }
  }
  // TODO: Here, I use the cell volumes from the neutronics model.  Instead, should
  // I sum the volumes of the local elements in step 3?
  for (const auto& kv : T_dot_V) {
    neutronics.set_temperature(kv.first, kv.second / neutronics.get_volume(kv.first));
  }
}

void CoupledDriver::update_density(bool relax)
{
  comm_.message("Updating density");
  auto& neutronics = this->get_neutronics_driver();
  auto& heat = this->get_heat_driver();

  // Step 1: On heat ranks, assign the current iterate of rho to the previous
  if (relax && heat.active()) {
    std::copy(cell_densities_.cbegin(), cell_densities_.cend(), cell_densities_.begin());
  }

  // Step 2: On heat ranks, get cell rho
  if (heat.active()) {
    auto elem_densities = heat.density_local();
    std::cout << "Rank, min dens, max dens: " << comm_.rank << ", "
              << *std::min_element(elem_densities.cbegin(), elem_densities.cend()) << ", "
              << *std::max_element(elem_densities.cbegin(), elem_densities.cend())
              << std::endl;
    heat.comm_.Barrier();

    for (gsl::index i = 0; i < cells_.size(); ++i) {
      if (cell_fluid_mask_.at(i) == 1) {
        double avg_rho = 0.0;
        double total_V = 0.0;
        for (const auto& e : cell_to_elems_.at(cells_.at(i))) {
          avg_rho += elem_densities.at(e) * elem_volumes_.at(e);
          total_V += elem_volumes_.at(e);
        }
        avg_rho /= total_V;
        Ensures(avg_rho > 0.0);
        cell_densities_.at(i) = avg_rho;
      }
    }
    if (relax) {
      if (alpha_rho_ == ROBBINS_MONRO) {
        int n = i_picard_ + 1;
        cell_densities_ = cell_densities_ / n + (1. - 1. / n) * cell_densities_prev_;
      } else {
        cell_densities_ =
          alpha_rho_ * cell_densities_ + (1.0 - alpha_rho_) * cell_densities_prev_;
      }
    }
  }

  // Step 3: On neutron ranks, accumulate cell rho from all heat ranks
  std::map<CellHandle, double> rho_dot_V;
  decltype(cells_) cells_recv;
  decltype(cell_volumes_) cell_volumes_recv;
  decltype(cell_densities_) cell_densities_recv;
  decltype(cell_fluid_mask_) cell_fluid_mask_recv;

  for (const auto& heat_rank : heat_ranks_) {
    comm_.send_and_recv(cells_recv, neutronics_root_, cells_, heat_rank);
    neutronics.comm_.broadcast(cells_recv);

    comm_.send_and_recv(cell_volumes_recv, neutronics_root_, cell_volumes_, heat_rank);
    neutronics.comm_.broadcast(cell_volumes_recv);

    comm_.send_and_recv(
      cell_densities_recv, neutronics_root_, cell_densities_, heat_rank);
    neutronics.comm_.broadcast(cell_densities_recv);

    comm_.send_and_recv(
      cell_fluid_mask_recv, neutronics_root_, cell_fluid_mask_, heat_rank);

    if (neutronics.active()) {
      for (gsl::index i = 0; i < cells_recv.size(); ++i) {
        if (cell_fluid_mask_recv.at(i) == 1) {
          auto rho = cell_densities_recv.at(i);
          auto V = cell_volumes_recv.at(i);
          rho_dot_V[cells_recv.at(i)] += rho * V;
        }
      }
    }
  }
  // TODO: Here, I use the cell volumes from the neutronics model.  Instead, should
  // I sum the volumes of the local elements in step 3?
  for (const auto& kv : rho_dot_V) {
    neutronics.set_density(kv.first, kv.second / neutronics.get_volume(kv.first));
  }
}

void CoupledDriver::init_mappings()
{
  comm_.message("Initializing mappings");
  const auto& heat = this->get_heat_driver();
  auto& neutronics = this->get_neutronics_driver();

  // Send and recv buffers
  std::vector<Position> centroids_send, centroids_recv;
  decltype(elem_to_cell_) elem_to_cell_send;

  for (const auto& heat_rank : heat_ranks_) {
    // Set mapping of local elem --> global cell handle
    if (comm_.rank == heat_rank) {
      centroids_send = heat.centroid_local();
    }
    this->comm_.send_and_recv(
      centroids_recv, neutronics_root_, centroids_send, heat_rank);
    if (comm_.rank == neutronics_root_) {
      elem_to_cell_send = neutronics.find(centroids_recv);
    }
    this->comm_.send_and_recv(
      elem_to_cell_, heat_rank, elem_to_cell_send, neutronics_root_);
    comm_.Barrier();
  }
  if (heat.active()) {
    // Set mapping of global cell handle -> local elem
    for (gsl::index e = 0; e < elem_to_cell_.size(); ++e) {
      auto c = elem_to_cell_[e];
      cell_to_elems_[c].push_back(e); // Use [ ] instead of at() to insert new item
    }
    // Get list of all global cell handles
    for (const auto& kv : cell_to_elems_) {
      cells_.push_back(kv.first);
    }
  }
  // Debug
  for (const auto& rank : heat_ranks_) {
    if (comm_.rank == rank) {
      std::cout << "Rank: " << comm_.rank << std::endl;
      std::cout << "Cells to elems: ";
      for (const auto& kv : cell_to_elems_) {
        std::cout << kv.first << " ";
      }
      std::cout << std::endl;
      std::cout << "Cells: ";
      for (const auto& c : cells_) {
        std::cout << c << " ";
      }
      std::cout << std::endl;
    }
    comm_.Barrier();
  }
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

  // Every heat rank keeps track of its own local cell T
  if (heat.active()) {
    auto sz = static_cast<unsigned long>(cells_.size());
    cell_temperatures_.resize({sz});
    cell_temperatures_prev_.resize({sz});
  }

  if (temperature_ic_ == Initial::neutronics) {
    // Send and recv buffers
    decltype(cells_) cells_recv;
    decltype(cell_temperatures_) cell_temperatures_send;
    // The neutronics root sends cell T to each heat rank
    for (const auto& heat_rank : heat_ranks_) {
      comm_.send_and_recv(cells_recv, neutronics_root_, cells_, heat_rank);
      if (comm_.rank == neutronics_root_) {
        const auto sz = static_cast<unsigned long>(cells_recv.size());
        cell_temperatures_send.resize({sz});
        for (gsl::index i = 0; i < sz; ++i) {
          cell_temperatures_send.at(i) = neutronics.get_temperature(cells_recv.at(i));
        }
      }
      comm_.send_and_recv(
        cell_temperatures_, heat_rank, cell_temperatures_send, neutronics_root_);
    }
  } else if (temperature_ic_ == Initial::heat) {
    // * This sets l_cel_temps_ on the the coupling_root, based on the
    //   temperatures received from the heat solver.
    // * We do not want to apply underrelaxation here (and at this point,
    //   there is no previous iterate of temperature, anyway).
    update_temperature(false);
  }

  // In both cases, only temperatures_ was set, so we explicitly set temperatures_prev_
  if (heat.active()) {
    std::copy(cell_temperatures_.begin(),
              cell_temperatures_.end(),
              cell_temperatures_prev_.begin());
  }
}

void CoupledDriver::init_volumes()
{
  comm_.message("Initializing volumes");
  const auto& heat = this->get_heat_driver();
  const auto& neutronics = this->get_neutronics_driver();

  if (heat.active()) {
    elem_volumes_ = heat.volume_local();
    for (const auto& c : cells_) {
      double V = 0.0;
      for (const auto& e : cell_to_elems_.at(c)) {
        V += elem_volumes_.at(e);
      }
      cell_volumes_.push_back(V);
    }
  }

#ifndef NDEBUG
  // Volume check.  2020-10-06: This is pretty expensive now, so it's only in DEBUG builds

  // An array of global cell volumes, which will be accumulated from local cell volumes.
  std::map<CellHandle, double> glob_volumes;

  // Get all local cell volumes from heat ranks and sum them into the global cell volumes.
  for (const auto& heat_rank : heat_ranks_) {
    decltype(cells_) cells_recv;
    decltype(cell_volumes_) cell_volumes_recv;
    comm_.send_and_recv(cells_recv, neutronics_root_, cells_, heat_rank);
    comm_.send_and_recv(cell_volumes_recv, neutronics_root_, cell_volumes_, heat_rank);
    if (comm_.rank == neutronics_root_) {
      for (gsl::index i = 0; i < cells_recv.size(); ++i) {
        glob_volumes[cells_recv.at(i)] += cell_volumes_recv.at(i);
      }
    }
  }
  comm_.Barrier();

  if (comm_.rank == neutronics_root_) {
    // Compare volume from neutron driver to accumulated volume
    for (const auto& kv : glob_volumes) {
      auto v_neutronics = neutronics.get_volume(kv.first);
      std::stringstream msg;
      msg << "Cell " << neutronics.cell_label(kv.first) << ", V = " << v_neutronics
          << " (Neutronics), " << kv.second << " (Accumulated from Heat/Fluids)";
      comm_.message(msg.str());
    }
  }
  comm_.Barrier();
#endif
}

void CoupledDriver::init_densities()
{
  comm_.message("Initializing densities");
  const auto& neutronics = this->get_neutronics_driver();
  const auto& heat = this->get_heat_driver();

  if (heat.active()) {
    auto sz = static_cast<unsigned long>(cells_.size());
    cell_densities_.resize({sz});
    cell_densities_prev_.resize({sz});
  }

  if (density_ic_ == Initial::neutronics) {
    decltype(cells_) cells_recv;
    decltype(cell_densities_) cell_densities_send;
    for (const auto& heat_rank : heat_ranks_) {
      comm_.send_and_recv(cells_recv, neutronics_root_, cells_, heat_rank);
      if (comm_.rank == neutronics_root_) {
        const auto sz = static_cast<unsigned long>(cells_recv.size());
        cell_densities_send.resize({sz});
        for (gsl::index i = 0; i < sz; ++i) {
          cell_densities_send.at(i) = neutronics.get_density(cells_recv.at(i));
        }
      }
      comm_.send_and_recv(
        cell_densities_, heat_rank, cell_densities_send, neutronics_root_);
    }
  } else if (density_ic_ == Initial::heat) {
    // * We do not want to apply underrelaxation here (and at this point,
    //   there is no previous iterate of density, anyway).
    update_density(false);
  }

  // In both cases, we need to explicitly set densities_prev_
  if (heat.active()) {
    std::copy(
      cell_densities_.cbegin(), cell_densities_.cend(), cell_densities_prev_.begin());
  }
}

void CoupledDriver::init_cell_fluid_mask()
{
  comm_.message("Initializing cell fluid mask");
  auto& heat = this->get_heat_driver();

  if (heat.active()) {
    auto elem_fluid_mask = heat.fluid_mask_local();
    for (const auto& kv : cell_to_elems_) {

      auto& elems = kv.second;
      auto in_fluid = elem_fluid_mask.at(elems.at(0));
      for (gsl::index i = 1; i < elems.size(); ++i) {
        if (in_fluid != elem_fluid_mask.at(elems.at(i))) {
          throw std::runtime_error("ENRICO detected a neutronics cell contains both "
                                   "fluid and solid T/H elements.");
        }
      }
      cell_fluid_mask_.push_back(in_fluid);
    }
  }
}

void CoupledDriver::init_heat_source()
{
  comm_.message("Initializing heat source");

  if (this->heat_fluids_driver_->active()) {
    auto sz = {cells_.size()};
    cell_heat_ = xt::empty<double>(sz);
    cell_heat_prev_ = xt::empty<double>(sz);
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
