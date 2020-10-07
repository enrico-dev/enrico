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

  // Discover the ranks that are in each comm
  neutronics_ranks_ = gather_subcomm_ranks(comm_, neutronics_comm);
  heat_ranks_ = gather_subcomm_ranks(comm_, neutronics_comm);

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
    throw std::runtime_error{"nek5000 was specified as a solver, but is not enabled in this build of ENRICO"};
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
      double l_sum = xt::eval(xt::sum(xt::abs(l_cell_temps_ - l_cell_temps_prev_)))[0];
      MPI_Reduce(&l_sum, &tnorm, 1, MPI_DOUBLE, MPI_SUM, 0, heat.comm_.comm);
      break;
    }
    case Norm::L2: {
      double l_sum = xt::eval(xt::sum(xt::square(l_cell_temps_ - l_cell_temps_prev_)))[0];
      double g_sum;
      MPI_Reduce(&l_sum, &g_sum, 1, MPI_DOUBLE, MPI_SUM, 0, heat.comm_.comm);
      tnorm = sqrt(g_sum);
      break;
    }
    case Norm::LINF: {
      double l_max = xt::eval(xt::amax(xt::abs(l_cell_temps_ - l_cell_temps_prev_)))[0];
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

void CoupledDriver::update_temperature(bool relax)
{
  comm_.message("Updating temperature");

  auto& neutronics = this->get_neutronics_driver();
  auto& heat = this->get_heat_driver();

  // Step 1: On each heat rank, assign the current iterate of T to the previous
  if (relax && heat.active()) {
    std::copy(l_cell_temps_.begin(), l_cell_temps_.end(), l_cell_temps_prev_.begin());
  }

  // Step 2: On each heat rank, compute the local cell T by volume-averaging the local
  // elem T
  if (heat.active()) {
    auto l_elem_temps = heat.temperature_local();
    for (gsl::index l_cell = 0; l_cell < n_local_cells_; ++l_cell) {
      double avg_temp = 0.0;
      double total_vol = 0.0;
      for (const auto& l_elem : l_cell_to_l_elems.at(l_cell)) {
        double T = l_elem_temps.at(l_elem);
        double V = l_elem_volumes_.at(l_elem);
        avg_temp += T * V;
        total_vol += V;
      }
      avg_temp /= total_vol;
      Ensures(avg_temp > 0.0);
      l_cell_temps_.at(l_cell) = avg_temp;
    }
    // Apply relaxation to local cell temperatures
    if (relax) {
      if (alpha_T_ == ROBBINS_MONRO) {
        int n = i_picard_ + 1;
        l_cell_temps_ = l_cell_temps_ / n + (1. - 1. / n) * l_cell_temps_prev_;
      } else {
        l_cell_temps_ =
          alpha_T_ * l_cell_temps_prev_ + (1.0 - alpha_T_) * l_cell_temps_prev_;
      }
    }
  }

  // Step 3: Set the global cell T to a null sentinel
  if (neutronics.active()) {
    for (CellHandle g_cell = 0; g_cell < neutronics.n_cells(); ++g_cell) {
      neutronics.set_temperature(g_cell, -1.0);
    }
  }

  // Step 4: Set the global cell T from the volume-averaged contributions the local cell T
  for (const auto& heat_rank : heat_ranks_) {
    comm_.send_and_recv(l_cell_temps_, neutronics_root_, heat_rank);
    neutronics.comm_.broadcast(l_cell_temps_);

    comm_.send_and_recv(l_cell_to_g_cell_, neutronics_root_, heat_rank);
    neutronics.comm_.broadcast(l_cell_to_g_cell_);

    comm_.send_and_recv(l_cell_volumes_, neutronics_root_, heat_rank);
    neutronics.comm_.broadcast(l_cell_volumes_);

    if (neutronics.active()) {
      for (gsl::index l_cell = 0; l_cell < l_cell_to_g_cell_.size(); ++l_cell) {
        auto g_cell = l_cell_to_g_cell_.at(l_cell);
        auto g_cell_T = neutronics.get_temperature(g_cell);
        auto l_cell_T = l_cell_temps_.at(l_cell);
        // If the global cell T is the null sentinel, then set the global cell
        // T to be the local cell T.
        if (g_cell_T < 0.0) {
          neutronics.set_temperature(g_cell, l_cell_T);
        }
        // Otherwise, add the contribution of the local cell (which was a part of the
        // global cell) to the existing temperature of the global cell
        else {
          auto g_cell_V = neutronics.get_volume(g_cell);
          auto l_cell_V = l_cell_volumes_.at(l_cell);
          auto avg_T =
            (g_cell_T * (g_cell_V - l_cell_V) + l_cell_T * l_cell_V) / g_cell_V;
          neutronics.set_temperature(g_cell, avg_T);
        }
      }
    }
  }

#ifndef NDEBUG
  // Make sure none of the temperatures remain null
  if (neutronics.active()) {
    for (CellHandle g_cell; g_cell < neutronics.n_cells(); ++g_cell) {
      Ensures(neutronics.get_temperature(g_cell) > 0.0);
    }
  }
#endif
}

void CoupledDriver::update_density(bool relax)
{

  auto& neutronics = this->get_neutronics_driver();
  auto& heat = this->get_heat_driver();

  // *********************************************************************
  // Gather densities on heat root and apply underrelaxation
  // *********************************************************************

  if (relax && comm_.rank == heat_root_) {
    std::copy(densities_.begin(), densities_.end(), densities_prev_.begin());
  }

  densities_ = heat.density();

  if (relax && comm_.rank == heat_root_) {
    if (alpha_rho_ == ROBBINS_MONRO) {
      int n = i_picard_ + 1;
      densities_ = densities_ / n + (1. - 1. / n) * densities_prev_;
    } else {
      densities_ = alpha_rho_ * densities_ + (1.0 - alpha_rho_) * densities_prev_;
    }
  }

  // ************************************************************************
  // Send underrelaxed density to all neutron ranks and set cell temperatures
  // ************************************************************************

  comm_.send_and_recv(densities_, neutronics_root_, heat_root_);
  neutronics.comm_.broadcast(densities_);

  if (neutronics.active()) {
    // For each neutronics cell in a fluid cell, volume average the densities
    // and set in driver
    for (const auto& kv : cell_to_elems_) {
      CellHandle cell = kv.first;
      if (cell_fluid_mask_[cell] == 1) {
        // Get volume-averaged density for this cell instance
        double average_density = 0.0;
        double total_vol = 0.0;
        for (int e : kv.second) {
          average_density += densities_[e] * elem_volumes_[e];
          total_vol += elem_volumes_[e];
        }
        // Set density for cell instance
        average_density /= total_vol;
        Ensures(average_density > 0.0);
        neutronics.set_density(cell, average_density);
      }
    }
  }
}

void CoupledDriver::init_mappings()
{
  comm_.message("Initializing mappings");

  const auto& heat = this->get_heat_driver();
  auto& neutronics = this->get_neutronics_driver();

  for (const auto& heat_rank : heat_ranks_) {
    // Get the local centroids for this heat/fluids rank
    std::vector<Position> l_cents;
    if (heat.comm_.rank == heat_rank) {
      l_cents = heat.centroid_local();
    }

    // (1) Send local elem centroids to neutron root.
    // (2) Determine mapping of local elems -> global cells.
    // (3) Send mapping back heat/fluids rank.
    this->comm_.send_and_recv(l_cents, neutronics_root_, heat_rank);
    if (neutronics.comm_.is_root()) {
      l_elem_to_g_cell_ = neutronics.find(l_cents);
    }
    this->comm_.send_and_recv(l_elem_to_g_cell_, heat_rank, neutronics_root_);
  }

  if (heat.active()) {
    // Establish a mapping of local cell IDs -> global cell IDs. This ends up with the
    // unique global IDs in order.
    l_cell_to_g_cell_.resize(l_elem_to_g_cell_.size());
    std::sort(l_cell_to_g_cell_.begin(), l_cell_to_g_cell_.end());
    auto new_end = std::unique(l_cell_to_g_cell_.begin(), l_cell_to_g_cell_.end());
    l_cell_to_g_cell_.erase(new_end, l_cell_to_g_cell_.end());

    n_local_cells_ = l_cell_to_g_cell_.size();

    // Determine mapping of local cells -> local elems
    l_cell_to_l_elems.resize(n_local_cells_);

    for (int32_t l_elem = 0; l_elem < l_elem_to_g_cell_.size(); ++l_elem) {
      auto g_cell = l_elem_to_g_cell_.at(l_elem);
      auto loc =
        std::lower_bound(l_cell_to_g_cell_.cbegin(), l_cell_to_g_cell_.cend(), g_cell);
      auto l_cell = std::distance(l_cell_to_g_cell_.cbegin(), loc);

      l_cell_to_l_elems.at(l_cell).push_back(l_elem);
    }
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
    l_cell_temps_.resize({static_cast<unsigned long>(n_local_cells_)});
    l_cell_temps_prev_.resize({static_cast<unsigned long>(n_local_cells_)});
  }

  if (temperature_ic_ == Initial::neutronics) {
    // The neutronics root sends cell T to each heat rank
    for (const auto& heat_rank : heat_ranks_) {
      comm_.send_and_recv(l_cell_to_g_cell_, neutronics_root_, heat_rank);
      if (comm_.rank == neutronics_root_) {
        for (gsl::index l_cell = 0; l_cell < l_cell_to_g_cell_.size(); ++l_cell) {
          l_cell_temps_.at(l_cell) =
            neutronics.get_temperature(l_cell_to_g_cell_.at(l_cell));
        }
      }
      comm_.send_and_recv(l_cell_temps_, heat_rank, neutronics_root_);
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
    std::copy(l_cell_temps_.begin(), l_cell_temps_.end(), l_cell_temps_prev_.begin());
  }
}

void CoupledDriver::init_volumes()
{
  comm_.message("Initializing volumes");

  const auto& heat = this->get_heat_driver();
  const auto& neutronics = this->get_neutronics_driver();

  // Each heat rank: (1) Gets its local elem volumes; (2) accumulates its local cell
  // volumes.  The local cells are ordered the same way as the keys in g_cell_to_l_elems_
  if (heat.active()) {
    l_cell_volumes_.resize(n_local_cells_);
    l_elem_volumes_ = heat.volume_local();
    for (const auto& kv : l_cell_to_l_elems) {
      double l_cell_vol = 0.;
      for (const auto& l_elem : kv.second) { // kv.second is l_elems
        l_cell_vol += l_elem_volumes_.at(l_elem);
      }
      l_cell_volumes_.at(kv.first) = l_cell_vol;
    }
  }

#ifndef NDEBUG
  // Volume check.  2020-10-06: This is pretty expensive now, so it's only in DEBUG builds

  // An array of global cell volumes, which will be accumulated from local cell volumes.
  std::vector<double> g_cell_vols;
  if (neutronics.comm_.is_root()) {
    g_cell_vols.resize(neutronics.n_cells());
  }

  // Get all local cell volumes from heat ranks and sum them into the global cell volumes.
  for (const auto& heat_rank : heat_ranks_) {
    comm_.send_and_recv(l_cell_to_g_cell_, neutronics_root_, heat_rank);
    comm_.send_and_recv(l_cell_volumes_, neutronics_root_, heat_rank);

    if (neutronics.comm_.is_root()) {
      assert(l_cell_to_g_cell_.size() == l_cell_volumes_.size());
      for (CellHandle l_cell = 0; l_cell < l_cell_to_g_cell_.size(); ++l_cell) {
        auto g_cell = l_cell_to_g_cell_.at(l_cell);
        g_cell_vols.at(g_cell) += l_cell_volumes_.at(l_cell);
      }
    }
  }

  // Compare volume from neutron driver to accumulated volume
  for (CellHandle g_cell; g_cell < neutronics.n_cells(); ++g_cell) {
    auto v_neutronics = neutronics.get_volume(g_cell);
    auto v_accum = g_cell_vols.at(g_cell);
    std::stringstream msg;
    msg << "Cell " << neutronics.cell_label(g_cell) << ", V = " << v_neutronics
        << " (Neutronics), " << v_accum << " (Accumulated from Heat/Fluids)";
    comm_.message(msg.str());
  }
#endif
}

void CoupledDriver::init_densities()
{
  comm_.message("Initializing densities");

  const auto& neutronics = this->get_neutronics_driver();
  const auto& heat = this->get_heat_driver();

  if (comm_.rank == heat_root_) {
    densities_.resize({static_cast<unsigned long>(n_global_elem_)});
    densities_prev_.resize({static_cast<unsigned long>(n_global_elem_)});
  } else if (comm_.rank == neutronics_root_) {
    densities_.resize({static_cast<unsigned long>(n_global_elem_)});
  }

  if (density_ic_ == Initial::neutronics) {
    if (comm_.rank == neutronics_root_) {
      // Loop over the neutronics cells, then loop over the TH elements
      // corresponding to that cell and assign the cell density to the correct
      // index in the densities_ array. This mapping assumes that each TH
      // element is fully contained within a neutronics cell, i.e., TH elements
      // are not split between multiple neutronics cells.
      for (const auto& kv : cell_to_elems_) {
        CellHandle cell = kv.first;
        const auto& global_elems = kv.second;

        for (int elem : global_elems) {
          if (cell_fluid_mask_[cell] == 1) {
            double rho = neutronics.get_density(cell);
            densities_[elem] = rho;
          } else {
            densities_[elem] = 0.0;
          }
        }
      }
    }
    comm_.send_and_recv(densities_, heat_root_, neutronics_root_);
  } else if (density_ic_ == Initial::heat) {
    // * This sets densities_ on the the coupling_root, based on the
    //   densities received from the heat solver.
    // * We do not want to apply underrelaxation here (and at this point,
    //   there is no previous iterate of density, anyway).
    update_density(false);
  }

  // In both cases, we need to explicitly set densities_prev_
  if (comm_.rank == heat_root_) {
    std::copy(densities_.begin(), densities_.end(), densities_prev_.begin());
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

  auto& heat = this->get_heat_driver();

  if (heat.active()) {
    l_cell_fluid_mask_.resize(n_local_cells_);
    auto l_elem_fluid_mask = heat.fluid_mask_local();

    for (const auto& kv : g_cell_to_l_elems_) {
      auto l_elem = kv.second.at(0);

      auto l_cell = g_cell_to_l_cell_.at(kv.first);
    }
  }

  for (gsl::index l_elem = 0; l_elem < l_elem_fluid_mask.size(); ++l_elem) {
  }

  // /* OLD ===================================================

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
    heat_source_ = xt::empty<double>({n_local_cells_});
    heat_source_prev_ = xt::empty<double>({n_local_cells_});
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
