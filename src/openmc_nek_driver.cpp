#include "enrico/openmc_nek_driver.h"

#include "enrico/const.h"
#include "enrico/error.h"
#include "enrico/message_passing.h"

#include "gsl/gsl"
#include "nek5000/core/nek_interface.h"
#include "openmc/capi.h"
#include "openmc/cell.h"
#include "pugixml.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xtensor.hpp"

#include <string>
#include <vector>

namespace enrico {

using gsl::index;

OpenmcNekDriver::OpenmcNekDriver(MPI_Comm comm, pugi::xml_node node)
  : CoupledDriver{comm, node}
{
  // Get parameters from enrico.xml
  pugi::xml_node nek_node = node.child("nek5000");
  double pressure_bc = node.child("pressure_bc").text().as_double();
  openmc_procs_per_node_ = node.child("openmc_procs_per_node").text().as_int();

  // Postcondition checks on user inputs
  Expects(openmc_procs_per_node_ > 0);

  // Create communicator for OpenMC with requested processes per node
  MPI_Comm openmc_comm;
  MPI_Comm intranode_comm;
  enrico::get_node_comms(
    comm_.comm, openmc_procs_per_node_, &openmc_comm, &intranode_comm);

  // Set intranode communicator
  intranode_comm_ = Comm(intranode_comm);

  // Instantiate OpenMC and Nek drivers
  openmc_driver_ = std::make_unique<OpenmcDriver>(openmc_comm);
  nek_driver_ = std::make_unique<NekDriver>(comm, pressure_bc, nek_node);

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

bool OpenmcNekDriver::has_global_coupling_data() const
{
  return openmc_driver_->active() && nek_driver_->active();
}

NeutronicsDriver& OpenmcNekDriver::get_neutronics_driver() const
{
  return *openmc_driver_;
}

HeatFluidsDriver& OpenmcNekDriver::get_heat_driver() const
{
  return *nek_driver_;
}

void OpenmcNekDriver::init_mappings()
{
  comm_.message("Initializing mappings");

  const auto& heat = this->get_heat_driver();
  if (heat.active()) {
    // Get centroids from heat driver
    auto elem_centroids = heat.centroids();

    // Broadcast centroids onto all the neutronics procs
    this->get_neutronics_driver().broadcast(elem_centroids);

    // Step 2: Set element->cell and cell->element mappings
    // Create buffer to store cell instance indices corresponding to each Nek global
    // element. This is needed because calls to OpenMC API functions can only be
    // made from processes
    std::vector<int32_t> elem_to_cell(heat.n_global_elem());

    if (openmc_driver_->active()) {
      std::unordered_map<CellInstance, index> cell_index;

      for (int32_t i = 0; i < elem_to_cell.size(); ++i) {
        // Determine cell instance corresponding to global element
        Position elem_pos = elem_centroids[i];
        CellInstance c{elem_pos};

        // If this cell instance hasn't been saved yet, add it to cells_ and
        // keep track of what index it corresponds to
        if (cell_index.find(c) == cell_index.end()) {
          cell_index[c] = openmc_driver_->cells_.size();
          openmc_driver_->cells_.push_back(c);
        }
        auto i_cell = cell_index.at(c);
        cell_to_elems_[i_cell].push_back(i);

        // Set value for cell instance in array
        elem_to_cell[i] = i_cell;
      }

      // Determine number of OpenMC cell instances
      n_cells_ = gsl::narrow<int32_t>(openmc_driver_->cells_.size());
    }

    // Set element -> cell instance mapping on each Nek rank
    intranode_comm_.Bcast(elem_to_cell.data(), elem_to_cell.size(), MPI_INT32_T);
    elem_to_cell_ = elem_to_cell;

    // Broadcast number of cell instances
    intranode_comm_.Bcast(&n_cells_, 1, MPI_INT32_T);
  }
}

void OpenmcNekDriver::init_tallies()
{
  using gsl::narrow_cast;

  comm_.message("Initializing tallies");

  if (openmc_driver_->active()) {
    // Build vector of material indices
    std::vector<openmc::CellInstance> instances;
    for (const auto& c : openmc_driver_->cells_) {
      instances.push_back(
        {narrow_cast<index>(c.index_), narrow_cast<index>(c.instance_)});
    }
    openmc_driver_->create_tallies(instances);
  }
}

void OpenmcNekDriver::init_temperatures()
{
  comm_.message("Initializing temperatures");

  if (this->has_global_coupling_data()) {
    auto n_global = this->get_heat_driver().n_global_elem();
    temperatures_.resize({n_global});
    temperatures_prev_.resize({n_global});

    if (temperature_ic_ == Initial::neutronics) {
      // Loop over the OpenMC cells, then loop over the global Nek elements
      // corresponding to that cell and assign the OpenMC cell temperature to
      // the correct index in the temperatures_ array. This mapping assumes that
      // each Nek element is fully contained within an OpenMC cell, i.e. Nek elements
      // are not split between multiple OpenMC cells.
      for (index i = 0; i < openmc_driver_->cells_.size(); ++i) {
        const auto& global_elems = cell_to_elems_.at(i);
        const auto& c = openmc_driver_->cells_[i];

        for (auto elem : global_elems) {
          double T = c.get_temperature();
          temperatures_[elem] = T;
          temperatures_prev_[elem] = T;
        }
      }
    } else if (temperature_ic_ == Initial::heat) {
      // Use whatever temperature is in Nek's internal arrays, either from a restart
      // file or from a useric fortran routine.
      update_temperature();

      // update_temperature() begins by saving temperatures_ to temperatures_prev_, and
      // then changes temperatures_. We need to save temperatures_ here to
      // temperatures_prev_ manually because init_temperatures() initializes both
      // temperatures_ and temperatures_prev_.
      std::copy(temperatures_.begin(), temperatures_.end(), temperatures_prev_.begin());
    }
  }
}

void OpenmcNekDriver::init_volumes()
{
  comm_.message("Initializing volumes");

  const auto& heat = this->get_heat_driver();
  if (heat.active()) {
    // Gather all the local element volumes on the Nek5000/OpenMC root
    elem_volumes_ = heat.volumes();

    // Broadcast global_element_volumes onto all the OpenMC procs
    this->get_neutronics_driver().broadcast(elem_volumes_);
  }

  // Volume check
  if (this->has_global_coupling_data()) {
    for (index i = 0; i < openmc_driver_->cells_.size(); ++i) {
      const auto& c = openmc_driver_->cells_[i];
      double v_openmc = c.volume_;
      double v_nek = 0.0;
      for (const auto& elem : cell_to_elems_.at(i)) {
        v_nek += elem_volumes_.at(elem);
      }
      std::stringstream msg;
      msg << "Cell " << openmc::model::cells[c.index_]->id_ << " (" << c.instance_
          << "), V = " << v_openmc << " (OpenMC), " << v_nek << " (Nek)";
      comm_.message(msg.str());
    }
  }
}

void OpenmcNekDriver::init_densities()
{
  comm_.message("Initializing densities");

  if (this->has_global_coupling_data()) {
    auto n_global = this->get_heat_driver().n_global_elem();
    densities_.resize({n_global});
    densities_prev_.resize({n_global});

    if (density_ic_ == Initial::neutronics) {
      // Loop over the OpenMC cells, then loop over the global Nek elements
      // corresponding to that cell and assign the OpenMC cell density to
      // the correct index in the densities_ array. This mapping assumes that
      // each Nek element is fully contained within an OpenMC cell, i.e. Nek
      // elements are not split between multiple OpenMC cells.
      for (index i = 0; i < openmc_driver_->cells_.size(); ++i) {
        auto& c = openmc_driver_->cells_[i];
        const auto& global_elems = cell_to_elems_.at(i);

        if (cell_fluid_mask_[i] == 1) {
          for (int elem : global_elems) {
            double rho = c.get_density();
            densities_[elem] = rho;
            densities_prev_[elem] = rho;
          }
        } else {
          for (int elem : global_elems) {
            densities_[elem] = 0.0;
            densities_prev_[elem] = 0.0;
          }
        }
      }
    } else if (density_ic_ == Initial::heat) {
      // Use whatever density is in Nek's internal arrays, either from a restart
      // file or from a useric fortran routine
      update_density();

      // update_density() begins by saving densities_ to densities_prev_, and
      // then changes densities_. We need to save densities_ here to densities_prev_
      // manually because init_densities() initializes both densities_ and
      // densities_prev_
      std::copy(densities_.begin(), densities_.end(), densities_prev_.begin());
    }
  }
}

void OpenmcNekDriver::init_elem_fluid_mask()
{
  comm_.message("Initializing element fluid mask");

  const auto& heat = this->get_heat_driver();

  if (heat.active()) {
    // On Nek's master rank, fm gets global data. On Nek's other ranks, fm is empty
    elem_fluid_mask_ = heat.fluid_mask();

    // Broadcast fluid mask to neutronics driver
    this->get_neutronics_driver().broadcast(elem_fluid_mask_);
  }
}

void OpenmcNekDriver::init_cell_fluid_mask()
{
  comm_.message("Initializing cell fluid mask");

  if (this->has_global_coupling_data()) {
    auto& cells = openmc_driver_->cells_;
    cell_fluid_mask_.resize({cells.size()});

    for (index i = 0; i < cells.size(); ++i) {
      auto elems = cell_to_elems_.at(i);
      for (const auto& j : elems) {
        if (elem_fluid_mask_[j] == 1) {
          cell_fluid_mask_[i] = 1;
          break;
        }
        cell_fluid_mask_[i] = 0;
      }
    }
  }
}

void OpenmcNekDriver::init_heat_source()
{
  comm_.message("Initializing heat source");

  heat_source_ = xt::empty<double>({n_cells_});
  heat_source_prev_ = xt::empty<double>({n_cells_});
}

void OpenmcNekDriver::set_heat_source()
{
  // OpenMC has heat source on each of its ranks. We need to make heat
  // source available on each Nek rank.
  intranode_comm_.Bcast(heat_source_.data(), n_cells_, MPI_DOUBLE);

  if (nek_driver_->active()) {
    // Determine displacement for this rank
    auto displacement = nek_driver_->local_displs_[nek_driver_->comm_.rank];

    // Loop over local elements to set heat source
    int n_local = this->get_heat_driver().n_local_elem();
    for (int32_t local_elem = 1; local_elem <= n_local; ++local_elem) {
      // get corresponding global element
      int32_t global_index = local_elem + displacement - 1;

      // get index to cell instance
      int32_t cell_index = elem_to_cell_.at(global_index);

      err_chk(nek_driver_->set_heat_source_at(local_elem, heat_source_[cell_index]),
              "Error setting heat source for local element " +
                std::to_string(local_elem));
    }
  }
}

void OpenmcNekDriver::set_temperature()
{
  if (this->get_heat_driver().active()) {
    if (openmc_driver_->active()) {
      // Broadcast global_element_temperatures onto all the OpenMC procs
      openmc_driver_->comm_.Bcast(temperatures_.data(), temperatures_.size(), MPI_DOUBLE);

      // For each OpenMC cell instance, volume average temperatures and set
      for (size_t i = 0; i < openmc_driver_->cells_.size(); ++i) {

        // Get corresponding global elements
        const auto& global_elems = cell_to_elems_.at(i);
        auto& c{openmc_driver_->cells_[i]};

        // Get volume-average temperature for this cell instance
        double average_temp = 0.0;
        double total_vol = 0.0;
        for (int elem : global_elems) {
          double T = temperatures_[elem];
          double V = elem_volumes_[elem];
          average_temp += T * V;
          total_vol += V;
        }

        // Set temperature for cell instance
        average_temp /= total_vol;
        Ensures(average_temp > 0.0);
        c.set_temperature(average_temp);
      }
    }
  }
}

void OpenmcNekDriver::set_density()
{
  if (this->get_heat_driver().active()) {
    // Since OpenMC's and Nek's master ranks are the same, we know that elem_densities_ on
    // OpenMC's master rank were updated.  Now we broadcast to the other OpenMC ranks.
    // TODO: This won't work if the Nek/OpenMC communicators are disjoint
    if (openmc_driver_->active()) {
      openmc_driver_->comm_.Bcast(densities_.data(), densities_.size(), MPI_DOUBLE);

      // For each OpenMC cell instance in a fluid cell, volume average the
      // densities and set
      // TODO:  Might be able to use xtensor masking to do some of this
      for (index i = 0; i < openmc_driver_->cells_.size(); ++i) {
        if (cell_fluid_mask_[i] == 1) {
          auto& c = openmc_driver_->cells_[i];
          double average_density = 0.0;
          double total_vol = 0.0;
          for (int e : cell_to_elems_.at(i)) {
            average_density += densities_[e] * elem_volumes_[e];
            total_vol += elem_volumes_[e];
          }
          double density = average_density / total_vol;
          Ensures(density > 0.0);
          c.set_density(average_density / total_vol);
        }
      }
    }
  }
}

} // namespace enrico
