#include "enrico/coupled_driver.h"

#include "enrico/comm_split.h"
#include "enrico/driver.h"
#include "enrico/error.h"
#include "enrico/nek_driver.h"
#include "enrico/openmc_driver.h"
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

  if (coup_node.child("alpha"))
    alpha_ = coup_node.child("alpha").text().as_double();

  if (coup_node.child("alpha_T"))
    alpha_T_ = coup_node.child("alpha_T").text().as_double();

  if (coup_node.child("alpha_rho"))
    alpha_rho_ = coup_node.child("alpha_rho").text().as_double();

  if (coup_node.child("temperature_ic")) {
    auto s = std::string{coup_node.child_value("temperature_ic")};

    if (s == "neutronics") {
      temperature_ic_ = Initial::neutronics;
    } else if (s == "heat_fluids") {
      temperature_ic_ = Initial::heat;
    } else {
      throw std::runtime_error{"Invalid value for <temperature_ic>"};
    }
  }

  if (coup_node.child("density_ic")) {
    auto s = std::string{coup_node.child_value("density_ic")};

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
  Expects(alpha_ > 0);
  Expects(alpha_T_ > 0);
  Expects(alpha_rho_ > 0);

  // Create communicators
  std::array<int, 2> nodes{neut_node.child("nodes").text().as_int(),
                           heat_node.child("nodes").text().as_int()};
  std::array<int, 2> procs_per_node{neut_node.child("procs_per_node").text().as_int(),
                                    heat_node.child("procs_per_node").text().as_int()};
  std::array<Comm, 2> driver_comms;
  Comm intranode_comm, coupling_comm; // Not used in current comm scheme

  get_driver_comms(
    comm_, nodes, procs_per_node, driver_comms, intranode_comm, coupling_comm);

  auto neutronics_comm = driver_comms[0];
  auto heat_comm = driver_comms[1];

  // Instantiate neutronics driver
  neutronics_driver_ = std::make_unique<OpenmcDriver>(neutronics_comm.comm);

  // Instantiate heat-fluids driver
  std::string s = heat_node.child_value("driver");
  if (s == "nek5000") {
    heat_fluids_driver_ = std::make_unique<NekDriver>(heat_comm.comm, heat_node);
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

  // Currently, the coupling root is the neutronics root.  Could be revisited
  coupling_root_ = neutronics_root_;

  // Send number of global elements to all procs
  n_global_elem_ = this->get_heat_driver().n_global_elem();
  comm_.broadcast(n_global_elem_, heat_root_);

  comm_.message("Communicator layout:");
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
      update_heat_source(i_timestep_ == 0 && i_picard_ == 0 ? false : true);

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

void CoupledDriver::compute_temperature_norm(const CoupledDriver::Norm& n,
                                             double& norm,
                                             bool& converged)
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

bool CoupledDriver::is_converged()
{
  bool converged;
  const auto& neutron = this->get_neutronics_driver();

  // assumes that the neutronics root has access to global coupling data
  if (comm_.rank == coupling_root_) {
    double norm;
    this->compute_temperature_norm(Norm::LINF, norm, converged);

    std::string msg = "temperature norm_linf: " + std::to_string(norm);
    neutron.comm_.message(msg);
  }

  comm_.broadcast(converged, coupling_root_);
  return converged;
}

void CoupledDriver::update_heat_source(bool relax)
{
  auto& neutronics = this->get_neutronics_driver();
  auto& heat = this->get_heat_driver();

  // ****************************************************************************
  // Send neutron driver's heat source to coupling root and apply underrelaxation
  // ****************************************************************************
  
  if (relax && comm_.rank == coupling_root_) {
      // Store previous heat source solution if more than one iteration has been performed
      // (otherwise there is not an initial condition for the heat source)
    std::copy(heat_source_.begin(), heat_source_.end(), heat_source_prev_.begin());
  }

  if (neutronics.active()) {
    heat_source_ = neutronics.heat_source(power_);
  }

  comm_.send_and_recv(heat_source_, coupling_root_, neutronics_root_);

  // Compute the next iterate of the heat source
  if (relax && comm_.rank == coupling_root_) {
    heat_source_ = alpha_ * heat_source_ + (1.0 - alpha_) * heat_source_prev_;
  }

  // ****************************************************************************
  // Send underrelaxed heat source to all heat ranks and set element temperatures
  // ****************************************************************************

  this->comm_.send_and_recv(heat_source_, heat_root_, coupling_root_);
  heat.comm_.broadcast(heat_source_);

  if (heat.active()) {
    // Determine displacement for this rank
    auto displacement = heat.local_displs_.at(heat.comm_.rank);
    int n_local_elem = heat.n_local_elem();
    // Set heat source in every element
    for (int32_t local_elem = 1; local_elem <= n_local_elem; ++local_elem) {
      int32_t global_elem = local_elem + displacement - 1;
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

  // *************************************************************************
  // Send heat driver's temperature to coupling root and apply underrelaxation
  // *************************************************************************
  
  if (relax && comm_.rank == coupling_root_) {
    // Store previous temperature solution; a previous solution will always be present
    // because a temperature IC is set and the neutronics solver runs first
    std::copy(temperatures_.begin(), temperatures_.end(), temperatures_prev_.begin());
  }

  temperatures_ = heat.temperature();
  this->comm_.send_and_recv(temperatures_, coupling_root_, heat_root_);

  if (relax && comm_.rank == coupling_root_) {
    temperatures_ = alpha_T_ * temperatures_ + (1.0 - alpha_T_) * temperatures_prev_;
  }

  // ****************************************************************************
  // Send underrelaxed temperature to all neutron ranks and set cell temperatures
  // ****************************************************************************

  // Broadcast global_element_temperatures onto all the neutronics procs
  comm_.send_and_recv(temperatures_, neutronics_root_, coupling_root_);
  neutronics.comm_.broadcast(temperatures_);

  if (neutronics.active()) {
    // For each neutronics cell, volume average temperatures and set
    for (CellHandle cell = 0; cell < cell_to_elems_.size(); ++cell) {
      // Get volume-average temperature for this cell instance
      double average_temp = 0.0;
      double total_vol = 0.0;
      for (int elem : cell_to_elems_.at(cell)) {
        double T = temperatures_[elem];
        double V = elem_volumes_[elem];
        average_temp += T * V;
        total_vol += V;
      }
      // Set temperature for cell instance
      average_temp /= total_vol;
      Ensures(average_temp > 0.0);
      neutronics.set_temperature(cell, average_temp);
    }
  }
}

void CoupledDriver::update_density(bool relax)
{

  auto& neutronics = this->get_neutronics_driver();
  auto& heat = this->get_heat_driver();

  // *********************************************************************
  // Send heat driver's density to coupling root and apply underrelaxation
  // *********************************************************************

  if (relax && comm_.rank == coupling_root_) {
    std::copy(densities_.begin(), densities_.end(), densities_prev_.begin());
  }

  densities_ = heat.density();
  comm_.send_and_recv(densities_, coupling_root_, heat_root_);

  if (relax && comm_.rank == coupling_root_) {
    densities_ = alpha_rho_ * densities_ (1.0 - alpha_rho_) * densities_prev_;
  }

  // ************************************************************************
  // Send underrelaxed density to all neutron ranks and set cell temperatures
  // ************************************************************************

  comm_.send_and_recv(densities_, neutronics_root_, coupling_root_);
  neutronics.comm_.broadcast(densities_);

  if (neutronics.active()) {
    // For each neutronics cell in a fluid cell, volume average the densities
    // and set in driver
    for (CellHandle cell = 0; cell < cell_to_elems_.size(); ++cell) {
      if (cell_fluid_mask_[cell] == 1) {
        // Get volume-averaged density for this cell instance
        double average_density = 0.0;
        double total_vol = 0.0;
        for (int e : cell_to_elems_.at(cell)) {
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

  std::vector<Position> elem_centroids;
  if (comm_.rank == heat_root_ || neutronics.active()) {
    elem_centroids.resize(n_global_elem_);
  }

  // Get centroids from heat driver and send to all neutronics procs
  elem_centroids = heat.centroids();
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

  if (comm_.rank == coupling_root_) {
    temperatures_.resize({static_cast<unsigned long>(n_global_elem_)});
    temperatures_prev_.resize({static_cast<unsigned long>(n_global_elem_)});
  }
  else if (comm_.rank == neutronics_root_) {
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
    comm_.send_and_recv(temperatures_, coupling_root_, neutronics_root_);
  } 
  else if (temperature_ic_ == Initial::heat) {
    // * This sets temperatures_ on the the coupling_root, based on the 
    //   temperatures received from the heat solver.
    // * We do not want to apply underrelaxation here (and at this point,
    //   there is no previous iterate of temperature, anyway).  
    update_temperature(false);
  }

  // In both cases, only temperatures_ was set, so we explicitly set temperatures_prev_
  if (comm_.rank == coupling_root_) {
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
    for (CellHandle cell = 0; cell < cell_to_elems_.size(); ++cell) {
      double v_neutronics = neutronics.get_volume(cell);
      double v_heatfluids = 0.0;
      for (const auto& elem : cell_to_elems_.at(cell)) {
        v_heatfluids += elem_volumes_.at(elem);
      }

      // TODO: Refactor to avoid dynamic_cast
      const auto* openmc_driver = dynamic_cast<const OpenmcDriver*>(&neutronics);
      if (openmc_driver) {
        const auto& c = openmc_driver->cells_[cell];
        std::stringstream msg;
        msg << "Cell " << openmc::model::cells[c.index_]->id_ << " (" << c.instance_
            << "), V = " << v_neutronics << " (Neutronics), " << v_heatfluids
            << " (Heat/Fluids)";
        comm_.message(msg.str());
      }
    }
  }
}

void CoupledDriver::init_densities()
{
  comm_.message("Initializing densities");

  const auto& neutronics = this->get_neutronics_driver();
  const auto& heat = this->get_heat_driver();

  if (comm_.rank == coupling_root_) {
    densities_.resize({static_cast<unsigned long>(n_global_elem_)});
    densities_prev_.resize({static_cast<unsigned long>(n_global_elem_)});
  } 
  else if (comm_.rank == neutronics_root_) {
    densities_.resize({static_cast<unsigned long>(n_global_elem_)});
  }

  if (density_ic_ == Initial::neutronics) {
    if (comm_.rank == neutronics_root_) {
      // Loop over the neutronics cells, then loop over the TH elements
      // corresponding to that cell and assign the cell density to the correct
      // index in the densities_ array. This mapping assumes that each TH
      // element is fully contained within a neutronics cell, i.e., TH elements
      // are not split between multiple neutronics cells.
      for (CellHandle cell = 0; cell < cell_to_elems_.size(); ++cell) {
        for (int elem : cell_to_elems_.at(cell)) {
          if (cell_fluid_mask_[cell] == 1) {
            double rho = neutronics.get_density(cell);
            densities_[elem] = rho;
          } else {
            densities_[elem] = 0.0;
          }
        }
      }
    }
    comm_.send_and_recv(densities_, coupling_root_, neutronics_root_);
  } 
  else if (density_ic_ == Initial::heat) {
    // * This sets densities_ on the the coupling_root, based on the 
    //   densities received from the heat solver.
    // * We do not want to apply underrelaxation here (and at this point,
    //   there is no previous iterate of density, anyway).  
    update_density(false);
  }

  // In both cases, we need to explicitly set densities_prev_
  std::copy(densities_.begin(), densities_.end(), densities_prev_.begin());
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

    for (CellHandle cell = 0; cell < n; ++cell) {
      auto elems = cell_to_elems_.at(cell);
      for (const auto& elem : elems) {
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

  if (comm_.rank == coupling_root_ || this->get_heat_driver().active()) {
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
