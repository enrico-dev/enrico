#include "enrico/coupled_driver.h"

#include "enrico/comm_split.h"
#include "enrico/driver.h"
#include "enrico/error.h"

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
#include <map>
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
  , timer_init_comms(comm_)
  , timer_init_mapping(comm_)
  , timer_init_tallies(comm_)
  , timer_init_volume(comm_)
  , timer_init_fluid_mask(comm_)
  , timer_init_temperature(comm_)
  , timer_init_density(comm_)
  , timer_init_heat_source(comm_)
  , timer_update_density(comm_)
  , timer_update_heat_source(comm_)
  , timer_update_temperature(comm_)
{
  parse_xml_params(node);
  init_comms(node);
  init_mapping();
  init_tallies();
  init_volume();
  init_fluid_mask();
  init_temperature();
  init_density();
  init_heat_source();
}

void CoupledDriver::parse_xml_params(const pugi::xml_node& node)
{
  auto coup_node = node.child("coupling");

  power_ = coup_node.child("power").text().as_double();
  max_timesteps_ = coup_node.child("max_timesteps").text().as_int();
  max_picard_iter_ = coup_node.child("max_picard_iter").text().as_int();
  verbose_ = coup_node.child("verbose").text().as_bool();
  if (coup_node.child("epsilon")) {
    epsilon_ = coup_node.child("epsilon").text().as_double();
  }

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
}

void CoupledDriver::init_comms(const pugi::xml_node& node)
{
  timer_init_comms.start();

  auto neut_node = node.child("neutronics");
  auto heat_node = node.child("heat_fluids");

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

  timer_init_comms.stop();

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

  timer_init_comms.start();

  // Discover the rank IDs (relative to comm_) that are in each single-physics subcomm
  neutronics_ranks_ = gather_subcomm_ranks(comm_, neutronics_comm);
  heat_ranks_ = gather_subcomm_ranks(comm_, heat_comm);

  // Send rank ID of neutronics subcomm root (relative to comm_) to all procs
  neutronics_root_ = this->get_neutronics_driver().comm_.is_root() ? comm_.rank : -1;
  MPI_Allreduce(MPI_IN_PLACE, &neutronics_root_, 1, MPI_INT, MPI_MAX, comm_.comm);

  // Send rank ID of heat subcomm root (relative to comm_) to all procs
  heat_root_ = this->get_heat_driver().comm_.is_root() ? comm_.rank : -1;
  MPI_Allreduce(MPI_IN_PLACE, &heat_root_, 1, MPI_INT, MPI_MAX, comm_.comm);

  timer_init_comms.stop();

  comm_report();
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
#ifdef _OPENMP
        omp_set_num_threads(neutronics.num_threads);
#pragma omp parallel default(none) shared(neutronics)
#pragma omp single
        {
          std::string msg = "OpenMP threads: " + std::to_string(omp_get_num_threads());
          neutronics.comm_.message(msg);
        }
#endif
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
#ifdef _OPENMP
        omp_set_num_threads(heat.num_threads);
#pragma omp parallel default(none) shared(heat)
#pragma omp single
        {
          std::string msg = "OpenMP threads: " + std::to_string(omp_get_num_threads());
          heat.comm_.message(msg);
        }
#endif
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

      timer_report();

      if (is_converged()) {
        std::string msg = "converged at i_picard = " + std::to_string(i_picard_);
        comm_.message(msg);
        break;
      }
    }
    comm_.Barrier();
  }
  // TODO: Is this final heat.write_step still needed?
  heat.write_step();
}

double CoupledDriver::temperature_norm(Norm norm)
{
  auto& heat = this->get_heat_driver();
  double global_norm = 0;

  // Update global temperature norm
  if (heat.active()) {
    switch (norm) {
    case Norm::L1: {
      double local_norm = xt::norm_l1(cell_temperature_ - cell_temperature_prev_)();
      MPI_Reduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, 0, heat.comm_.comm);
      break;
    }
    case Norm::L2: {
      double local_norm = xt::norm_sq(cell_temperature_ - cell_temperature_prev_)();
      MPI_Reduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, 0, heat.comm_.comm);
      global_norm = std::sqrt(global_norm);
      break;
    }
    case Norm::LINF: {
      double local_norm = xt::norm_linf(cell_temperature_ - cell_temperature_prev_)();
      MPI_Reduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_MAX, 0, heat.comm_.comm);
      break;
    }
    }
  }

  return global_norm;
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
  timer_update_heat_source.start();

  auto& neutronics = this->get_neutronics_driver();
  auto& heat = this->get_heat_driver();

  if (relax && heat.active()) {
    std::copy(
      cell_heat_source_.cbegin(), cell_heat_source_.cend(),
              cell_heat_source_prev_.begin());
  }

  decltype(cell_to_glob_cell_) cells_recv;
  decltype(cell_heat_source_) cell_heat_send;
  xt::xtensor<double, 1> all_cell_heat;

  // For the coupling scheme, only the neutronics root needs the heat source.
  // However, to compute the heat source, OpenmcDriver::heat_source must
  // do a collective operation on all the ranks in the neutronics sub comm.
  // Hence, all neutronics ranks must call OpenmcDriver::heat_source
  if (neutronics.active()) {
    all_cell_heat = neutronics.heat_source(power_);
  }

  // The neutronics root sends the cell-averaged heat sources to the heat ranks.
  // Each heat rank gets only the heat sources for its local cells.
  for (const auto& heat_rank : heat_ranks_) {
    comm_.send_and_recv(cells_recv, neutronics_root_, cell_to_glob_cell_, heat_rank);
    cell_heat_send.resize({cells_recv.size()});
    if (comm_.rank == neutronics_root_) {
      for (gsl::index i = 0; i < cells_recv.size(); ++i) {
        auto j = neutronics.cell_index(cells_recv.at(i));
        cell_heat_send.at(i) = all_cell_heat.at(j);
      }
    }
    comm_.send_and_recv(cell_heat_source_, heat_rank, cell_heat_send, neutronics_root_);
  }

  // On heat rank, update the elements' heat sources based on the cell-avged heat sources
  if (heat.active()) {
    if (relax) {
      if (alpha_ == ROBBINS_MONRO) {
        int n = i_picard_ + 1;
        cell_heat_source_ = cell_heat_source_ / n + (1. - 1. / n) * cell_heat_source_prev_;
      } else {
        cell_heat_source_ = alpha_ * cell_heat_source_ + (1.0 - alpha_) * cell_heat_source_prev_;
      }
    }
    for (gsl::index i = 0; i < cell_to_glob_cell_.size(); ++i) {
      for (const auto& e : glob_cell_to_elem_.at(cell_to_glob_cell_.at(i))) {
        heat.set_heat_source_at(e, cell_heat_source_.at(i));
      }
    }
  }
  timer_update_heat_source.stop();
}

void CoupledDriver::update_temperature(bool relax)
{
  comm_.message("Updating temperature");
  timer_update_temperature.start();

  auto& neutronics = this->get_neutronics_driver();
  auto& heat = this->get_heat_driver();

  // Step 1: On each heat rank, assign the current iterate of local cell-avged T
  // to the previous iterate
  if (relax && heat.active()) {
    std::copy(cell_temperature_.begin(),
              cell_temperature_.end(), cell_temperature_prev_.begin());
  }

  // Step 2: On each heat, compute cell-avged T
  if (heat.active()) {
    auto elem_temperatures = heat.temperature();
    for (gsl::index i = 0; i < cell_to_glob_cell_.size(); ++i) {
      double T_avg = 0.0;
      double V_tot = 0.0;
      for (const auto& e : glob_cell_to_elem_.at(cell_to_glob_cell_.at(i))) {
        double T = elem_temperatures.at(e);
        double V = elem_volume_.at(e);
        T_avg += T * V;
      }
      T_avg /= cell_volume_.at(i);
      Ensures(T_avg > 0.0);
      cell_temperature_.at(i) = T_avg;
    }
    // Apply relaxation to local cell-avged T
    if (relax) {
      if (alpha_T_ == ROBBINS_MONRO) {
        int n = i_picard_ + 1;
        cell_temperature_ =
          cell_temperature_ / n + (1. - 1. / n) * cell_temperature_prev_;
      } else {
        cell_temperature_ =
          alpha_T_ * cell_temperature_ + (1.0 - alpha_T_) * cell_temperature_prev_;
      }
    }
  }

  // Step 3: On each neutron rank, accumulate cell-avged volumes from all heat ranks
  std::unordered_map<CellHandle, double> T_dot_V;
  std::unordered_map<CellHandle, double> cell_V;
  decltype(cell_to_glob_cell_) cells_recv;
  decltype(cell_volume_) cell_volumes_recv;
  decltype(cell_temperature_) cell_temperatures_recv;
  for (const auto& heat_rank : heat_ranks_) {
    comm_.send_and_recv(cells_recv, neutronics_root_, cell_to_glob_cell_, heat_rank);
    neutronics.comm_.broadcast(cells_recv);

    comm_.send_and_recv(
      cell_temperatures_recv, neutronics_root_, cell_temperature_, heat_rank);
    neutronics.comm_.broadcast(cell_temperatures_recv);

    comm_.send_and_recv(cell_volumes_recv, neutronics_root_, cell_volume_, heat_rank);
    neutronics.comm_.broadcast(cell_volumes_recv);

    if (neutronics.active()) {
      for (gsl::index i = 0; i < cells_recv.size(); ++i) {
        auto T = cell_temperatures_recv.at(i);
        auto V = cell_volumes_recv.at(i);
        cell_V[cells_recv.at(i)] += V;
        T_dot_V[cells_recv.at(i)] += T * V;
      }
    }
  }
  for (const auto& kv : T_dot_V) {
    auto cell = kv.first;
    auto tv = kv.second;
    neutronics.set_temperature(cell, tv / cell_V.at(cell));
  }
  timer_update_temperature.stop();
}

void CoupledDriver::update_density(bool relax)
{
  comm_.message("Updating density");
  timer_update_density.start();

  auto& neutronics = this->get_neutronics_driver();
  auto& heat = this->get_heat_driver();

  // Step 1: On each heat rank, assign the current iterate of local cell-avged rho
  // to the previous iterate
  if (relax && heat.active()) {
    std::copy(cell_density_.cbegin(), cell_density_.cend(), cell_density_.begin());
  }

  // Step 2: On each heat, compute cell-avged rho
  if (heat.active()) {
    auto elem_densities = heat.density();

    for (gsl::index i = 0; i < cell_to_glob_cell_.size(); ++i) {
      if (cell_fluid_mask_.at(i) == 1) {
        double rho_avg = 0.0;
        double V_tot = 0.0;
        for (const auto& e : glob_cell_to_elem_.at(cell_to_glob_cell_.at(i))) {
          rho_avg += elem_densities.at(e) * elem_volume_.at(e);
        }
        rho_avg /= cell_volume_[i];
        Ensures(rho_avg > 0.0);
        cell_density_.at(i) = rho_avg;
      }
    }
    if (relax) {
      if (alpha_rho_ == ROBBINS_MONRO) {
        int n = i_picard_ + 1;
        cell_density_ = cell_density_ / n + (1. - 1. / n) * cell_density_prev_;
      } else {
        cell_density_ =
          alpha_rho_ * cell_density_ + (1.0 - alpha_rho_) * cell_density_prev_;
      }
    }
  }

  // Step 3: On each neutron rank, accumulate cell-avged volumes from all heat ranks
  std::map<CellHandle, double> rho_dot_V;
  std::map<CellHandle, double> cell_V;
  decltype(cell_to_glob_cell_) cells_recv;
  decltype(cell_volume_) cell_volumes_recv;
  decltype(cell_density_) cell_densities_recv;
  decltype(cell_fluid_mask_) cell_fluid_mask_recv;

  for (const auto& heat_rank : heat_ranks_) {
    comm_.send_and_recv(cells_recv, neutronics_root_, cell_to_glob_cell_, heat_rank);
    neutronics.comm_.broadcast(cells_recv);

    comm_.send_and_recv(cell_volumes_recv, neutronics_root_, cell_volume_, heat_rank);
    neutronics.comm_.broadcast(cell_volumes_recv);

    comm_.send_and_recv(
      cell_densities_recv, neutronics_root_, cell_density_, heat_rank);
    neutronics.comm_.broadcast(cell_densities_recv);

    comm_.send_and_recv(
      cell_fluid_mask_recv, neutronics_root_, cell_fluid_mask_, heat_rank);
    neutronics.comm_.broadcast(cell_fluid_mask_recv);

    if (neutronics.active()) {
      for (gsl::index i = 0; i < cells_recv.size(); ++i) {
        if (cell_fluid_mask_recv.at(i) == 1) {
          auto rho = cell_densities_recv.at(i);
          auto V = cell_volumes_recv.at(i);
          cell_V[cells_recv.at(i)] += V;
          rho_dot_V[cells_recv.at(i)] += rho * V;
        }
      }
    }
  }

  for (const auto& kv : rho_dot_V) {
    neutronics.set_density(kv.first, kv.second / cell_V.at(kv.first));
  }
  timer_update_density.stop();
}

void CoupledDriver::init_mapping()
{
  comm_.message("Initializing mappings");
  timer_init_mapping.start();

  const auto& heat = this->get_heat_driver();
  auto& neutronics = this->get_neutronics_driver();

  // Send and recv buffers
  std::vector<Position> centroids_send;
  std::vector<Position> centroids_recv;
  decltype(elem_to_glob_cell_) elem_to_cell_send;

  for (const auto& heat_rank : heat_ranks_) {
    // For the given heat rank, the neutronics root discovers the mapping of
    // local elem ID --> global cell handle.
    // * IMPORTANT: OpenmcDriver::find adds the cell instances it discovers to
    //   the OpenmcDriver::cells_ array of the calling rank only.  However, every
    //   neutronics rank needs the full array of cells_ for Openmc::create_tallies.
    //   Hence, we broadcast the centroids to all neutronics ranks and then
    //   call neutronics.find on each neutronics rank.
    if (comm_.rank == heat_rank) {
      centroids_send = heat.centroid();
    }
    this->comm_.send_and_recv(
      centroids_recv, neutronics_root_, centroids_send, heat_rank);
    neutronics.comm_.broadcast(centroids_recv);
    if (neutronics.comm_.active()) {
      elem_to_cell_send = neutronics.find(centroids_recv);
    }

    // The neutronics root send the mapping of local elem ID --> global cell handle
    // back to the given heat rank.
    this->comm_.send_and_recv(
      elem_to_glob_cell_, heat_rank, elem_to_cell_send, neutronics_root_);
    comm_.Barrier();
  }
  if (heat.active()) {
    // The heat rank sets the inverse mapping of global cell handle -> local element ID
    // This is only for its local cells.
    for (gsl::index e = 0; e < elem_to_glob_cell_.size(); ++e) {
      auto c = elem_to_glob_cell_[e];
      glob_cell_to_elem_[c].push_back(e); // Use [ ] instead of at() to insert new item
    }
    // The heat rank creates an array global cell handles for its local cells.
    // This is useful in the coupling.
    for (const auto& kv : glob_cell_to_elem_) {
      cell_to_glob_cell_.push_back(kv.first);
    }
  }
  timer_init_mapping.stop();
}

void CoupledDriver::init_tallies()
{
  comm_.message("Initializing tallies");
  timer_init_tallies.start();

  auto& neutronics = this->get_neutronics_driver();
  if (neutronics.active()) {
    neutronics.create_tallies();
  }
  timer_init_tallies.stop();
}

void CoupledDriver::init_temperature()
{
  comm_.message("Initializing temperatures");
  timer_init_temperature.start();

  const auto& neutronics = this->get_neutronics_driver();
  const auto& heat = this->get_heat_driver();

  // Every heat rank keeps track of its own local cell T
  if (heat.active()) {
    auto sz = static_cast<unsigned long>(cell_to_glob_cell_.size());
    cell_temperature_.resize({sz});
    cell_temperature_prev_.resize({sz});
  }

  if (temperature_ic_ == Initial::neutronics) {
    // Send and recv buffers
    decltype(cell_to_glob_cell_) cells_recv;
    decltype(cell_temperature_) cell_temperatures_send;
    // The neutronics root sends cell T to each heat rank
    for (const auto& heat_rank : heat_ranks_) {
      comm_.send_and_recv(cells_recv, neutronics_root_, cell_to_glob_cell_, heat_rank);
      if (comm_.rank == neutronics_root_) {
        const auto sz = static_cast<unsigned long>(cells_recv.size());
        cell_temperatures_send.resize({sz});
        for (gsl::index i = 0; i < sz; ++i) {
          cell_temperatures_send.at(i) = neutronics.get_temperature(cells_recv.at(i));
        }
      }
      comm_.send_and_recv(
        cell_temperature_, heat_rank, cell_temperatures_send, neutronics_root_);
    }
  } else if (temperature_ic_ == Initial::heat) {
    //  We do not want to apply underrelaxation here since, at this point, there is no
    //  previous iterate of temperature.
    update_temperature(false);
  }

  // In both cases, only temperatures_ was set, so we explicitly set temperatures_prev_
  if (heat.active()) {
    std::copy(cell_temperature_.begin(),
              cell_temperature_.end(), cell_temperature_prev_.begin());
  }
  timer_init_temperature.stop();
}

void CoupledDriver::init_volume()
{
  comm_.message("Initializing volumes");
  timer_init_volume.start();

  const auto& heat = this->get_heat_driver();
  const auto& neutronics = this->get_neutronics_driver();

  if (heat.active()) {
    elem_volume_ = heat.volume();
    for (const auto& c : cell_to_glob_cell_) {
      double V = 0.0;
      for (const auto& e : glob_cell_to_elem_.at(c)) {
        V += elem_volume_.at(e);
      }
      cell_volume_.push_back(V);
    }
  }
  timer_init_volume.stop();

  check_volumes();
}

void CoupledDriver::check_volumes()
{
  comm_.message("Volume check");
  const auto& neutronics = this->get_neutronics_driver();

  // An array of global cell volumes, which will be accumulated from local cell volumes.
  std::map<CellHandle, double> glob_volumes;

  // Get all local cell volumes from heat ranks and sum them into the global cell volumes.
  for (const auto& heat_rank : heat_ranks_) {
    decltype(cell_to_glob_cell_) cells_recv;
    decltype(cell_volume_) cell_volumes_recv;
    comm_.send_and_recv(cells_recv, neutronics_root_, cell_to_glob_cell_, heat_rank);
    comm_.send_and_recv(cell_volumes_recv, neutronics_root_, cell_volume_, heat_rank);
    if (comm_.rank == neutronics_root_) {
      for (gsl::index i = 0; i < cells_recv.size(); ++i) {
        glob_volumes[cells_recv.at(i)] += cell_volumes_recv.at(i);
      }
    }
  }
  comm_.Barrier();

  if (comm_.rank == neutronics_root_) {
    double v_rel_diff_min = std::numeric_limits<double>::max();
    double v_rel_diff_max = -1;
    double v_rel_diff_sum = 0;
    double v_rel_diff_count = 0;

    // Compare volume from neutron driver to accumulated volume
    for (const auto& kv : glob_volumes) {
      auto cell = kv.first;
      auto v_accum = kv.second;
      auto v_neutronics = neutronics.get_volume(cell);

      // In neutronics model, volume = 1.0 is a dummy value
      if (v_neutronics != 1.0) {
        auto v_diff = std::abs(v_neutronics - v_accum);
        auto v_rel_diff = v_diff / v_neutronics;

        v_rel_diff_min = std::min(v_rel_diff_min, v_rel_diff);
        v_rel_diff_max = std::max(v_rel_diff_max, v_rel_diff);
        v_rel_diff_sum += v_rel_diff;
        v_rel_diff_count += 1;

        if (verbose_) {
          std::stringstream msg;
          msg << "  Cell " << neutronics.cell_label(cell) << ", V = " << v_neutronics
              << " (Neutronics), " << v_accum << " (Accumulated from Heat/Fluids)";
          comm_.message(msg.str(), neutronics_root_);
        }
      }
    }

    std::stringstream msg;
    msg << "  Min relative volume diff:  " << v_rel_diff_min;
    comm_.message(msg.str(), neutronics_root_);

    msg.str("");
    msg << "  Max relative volume diff:  " << v_rel_diff_max;
    comm_.message(msg.str(), neutronics_root_);

    msg.str("");
    msg << "  Mean relative volume diff: " << v_rel_diff_sum / v_rel_diff_count;
    comm_.message(msg.str(), neutronics_root_);
  }
  comm_.Barrier();
}

void CoupledDriver::init_density()
{
  comm_.message("Initializing densities");
  timer_init_density.start();

  const auto& neutronics = this->get_neutronics_driver();
  const auto& heat = this->get_heat_driver();

  if (heat.active()) {
    auto sz = static_cast<unsigned long>(cell_to_glob_cell_.size());
    cell_density_.resize({sz});
    cell_density_prev_.resize({sz});
  }

  if (density_ic_ == Initial::neutronics) {
    decltype(cell_to_glob_cell_) cells_recv;
    decltype(cell_density_) cell_densities_send;
    for (const auto& heat_rank : heat_ranks_) {
      comm_.send_and_recv(cells_recv, neutronics_root_, cell_to_glob_cell_, heat_rank);
      if (comm_.rank == neutronics_root_) {
        const auto sz = static_cast<unsigned long>(cells_recv.size());
        cell_densities_send.resize({sz});
        for (gsl::index i = 0; i < sz; ++i) {
          cell_densities_send.at(i) = neutronics.get_density(cells_recv.at(i));
        }
      }
      comm_.send_and_recv(
        cell_density_, heat_rank, cell_densities_send, neutronics_root_);
    }
  } else if (density_ic_ == Initial::heat) {
    // * We do not want to apply underrelaxation here (and at this point,
    //   there is no previous iterate of density, anyway).
    update_density(false);
  }

  // In both cases, we need to explicitly set densities_prev_
  if (heat.active()) {
    std::copy(cell_density_.cbegin(), cell_density_.cend(), cell_density_prev_.begin());
  }
  timer_init_density.stop();
}

void CoupledDriver::init_fluid_mask()
{
  comm_.message("Initializing cell fluid mask");
  timer_init_fluid_mask.start();

  auto& heat = this->get_heat_driver();

  if (heat.active()) {
    auto elem_fluid_mask = heat.fluid_mask();
    for (const auto& kv : glob_cell_to_elem_) {

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
  timer_init_fluid_mask.stop();
}

void CoupledDriver::init_heat_source()
{
  comm_.message("Initializing heat source");
  timer_init_heat_source.start();

  if (this->heat_fluids_driver_->active()) {
    auto sz = {cell_to_glob_cell_.size()};
    cell_heat_source_ = xt::empty<double>(sz);
    cell_heat_source_prev_ = xt::empty<double>(sz);
  }
  timer_init_heat_source.stop();
}

void CoupledDriver::comm_report()
{
  if (!verbose_)
    return;

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

void CoupledDriver::timer_report()
{

  auto& heat = this->get_heat_driver();
  auto& neut = this->get_neutronics_driver();

  std::vector<TimeAmt> coup_times{
    {"init_comms", timer_init_comms.elapsed()},
    {"init_fluid_mask", timer_init_fluid_mask.elapsed()},
    {"init_density", timer_init_density.elapsed()},
    {"init_heat_source", timer_init_heat_source.elapsed()},
    {"init_mapping", timer_init_mapping.elapsed()},
    {"init_tallies", timer_init_tallies.elapsed()},
    {"init_temperature", timer_init_temperature.elapsed()},
    {"init_volume", timer_init_volume.elapsed()},
    {"update_density", timer_update_density.elapsed()},
    {"update_heat_source", timer_update_heat_source.elapsed()},
    {"update_temperature", timer_update_temperature.elapsed()}};

  std::vector<TimeAmt> heat_times{{"driver_setup", heat.timer_driver_setup.elapsed()},
                                  {"init_step", heat.timer_init_step.elapsed()},
                                  {"solve_step", heat.timer_solve_step.elapsed()},
                                  {"write_step", heat.timer_write_step.elapsed()},
                                  {"finalize_step", heat.timer_finalize_step.elapsed()}};

  std::vector<TimeAmt> neut_times{{"driver_setup", neut.timer_driver_setup.elapsed()},
                                  {"init_step", neut.timer_init_step.elapsed()},
                                  {"solve_step", neut.timer_solve_step.elapsed()},
                                  {"write_step", neut.timer_write_step.elapsed()},
                                  {"finalize_step", neut.timer_finalize_step.elapsed()}};

  auto tot_time = TimeAmt::sum_times(coup_times) + TimeAmt::sum_times(heat_times) +
                  TimeAmt::sum_times(neut_times);
  auto nrm = [tot_time](TimeAmt& t) { t.percent = t.time / tot_time; };
  std::for_each(coup_times.begin(), coup_times.end(), nrm);
  std::for_each(heat_times.begin(), heat_times.end(), nrm);
  std::for_each(neut_times.begin(), neut_times.end(), nrm);

  std::stringstream msg;
  msg << "Cumulative times at i_timestep = " << i_timestep_ << " , i_picard = " << i_picard_;
  comm_.message(msg.str());

  TimeAmt::print_times("CoupledDriver", coup_times, comm_);
  TimeAmt::print_times("NeutronicsDriver", neut_times, comm_);
  TimeAmt::print_times("HeatFluidsDriver", heat_times, comm_);

  auto tot_pct = TimeAmt::sum_percent(coup_times) + TimeAmt::sum_percent(heat_times) +
                 TimeAmt::sum_percent(neut_times);
  std::vector<TimeAmt> total_time{{"total", tot_time, tot_pct}};
  TimeAmt::print_times("Total", total_time, comm_);
}

} // namespace enrico
