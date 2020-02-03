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

    // Set element->cell and cell->element mappings. Create buffer to store cell
    // handles corresponding to each heat-fluids global element.
    std::vector<CellHandle> elem_to_cell(heat.n_global_elem());

    auto& neutronics = this->get_neutronics_driver();
    if (neutronics.active()) {
      // Get cell handle corresponding to each element centroid
      elem_to_cell = neutronics.find(elem_centroids);

      // Create a vector of elements for each neutronics cell
      for (int32_t elem = 0; elem < elem_to_cell.size(); ++elem) {
        auto cell = elem_to_cell[elem];
        cell_to_elems_[cell].push_back(elem);
      }

      // Determine number of OpenMC cell instances
      n_cells_ = cell_to_elems_.size();
    }

    // Set element -> cell instance mapping on each Nek rank
    intranode_comm_.Bcast(elem_to_cell.data(), elem_to_cell.size(), MPI_UINT64_T);
    elem_to_cell_ = elem_to_cell;

    // Broadcast number of cell instances
    intranode_comm_.Bcast(&n_cells_, 1, MPI_INT32_T);
  }
}

void OpenmcNekDriver::init_tallies()
{

  comm_.message("Initializing tallies");

  auto& neutronics = this->get_neutronics_driver();
  if (neutronics.active()) {
    auto n = cell_to_elems_.size();
    neutronics.create_tallies(n);
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
      // Loop over heat-fluids elements and determine temperature based on
      // corresponding neutronics cell. This mapping assumes that each
      // heat-fluids element is fully contained within a neutronic cell, i.e.,
      // heat-fluids elements are not split between multiple neutronics cells.
      const auto& neutronics = this->get_neutronics_driver();
      for (index elem = 0; elem < elem_to_cell_.size(); ++elem) {
        auto cell = elem_to_cell_[elem];
        double T = neutronics.get_temperature(cell);
        temperatures_[elem] = T;
        temperatures_prev_[elem] = T;
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
    const auto& neutronics = this->get_neutronics_driver();
    for (CellHandle cell = 0; cell < cell_to_elems_.size(); ++cell) {
      double v_openmc = neutronics.get_volume(cell);
      double v_nek = 0.0;
      for (const auto& elem : cell_to_elems_.at(cell)) {
        v_nek += elem_volumes_.at(elem);
      }

      // TODO: Refactor to avoid dynamic_cast
      const auto* openmc_driver = dynamic_cast<const OpenmcDriver*>(&neutronics);
      if (openmc_driver) {
        const auto& c = openmc_driver->cells_[cell];
        std::stringstream msg;
        msg << "Cell " << openmc::model::cells[c.index_]->id_ << " (" << c.instance_
            << "), V = " << v_openmc << " (OpenMC), " << v_nek << " (Nek)";
        comm_.message(msg.str());
      }
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
      const auto& neutronics = this->get_neutronics_driver();
      for (CellHandle cell = 0; cell < cell_to_elems_.size(); ++cell) {
        const auto& global_elems = cell_to_elems_.at(cell);

        if (cell_fluid_mask_[cell] == 1) {
          for (int elem : global_elems) {
            double rho = neutronics.get_density(cell);
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

  auto& heat = this->get_heat_driver();
  if (heat.active()) {
    // Determine displacement for this rank
    auto displacement = heat.local_displs_[heat.comm_.rank];

    // Loop over local elements to set heat source
    int n_local = heat.n_local_elem();
    for (int32_t local_elem = 1; local_elem <= n_local; ++local_elem) {
      // get corresponding global element
      int32_t global_elem = local_elem + displacement - 1;

      // get index to cell instance
      CellHandle cell = elem_to_cell_.at(global_elem);

      err_chk(heat.set_heat_source_at(local_elem, heat_source_[cell]),
              "Error setting heat source for local element " +
                std::to_string(local_elem));
    }
  }
}

void OpenmcNekDriver::set_temperature()
{
  if (this->get_heat_driver().active()) {
    auto& neutronics = this->get_neutronics_driver();
    if (neutronics.active()) {
      // Broadcast global_element_temperatures onto all the OpenMC procs
      neutronics.comm_.Bcast(temperatures_.data(), temperatures_.size(), MPI_DOUBLE);

      // For each OpenMC cell instance, volume average temperatures and set
      for (CellHandle cell = 0; cell < cell_to_elems_.size(); ++cell) {

        // Get corresponding global elements
        const auto& global_elems = cell_to_elems_.at(cell);

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
        neutronics.set_temperature(cell, average_temp);
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
    auto& neutronics = this->get_neutronics_driver();
    if (neutronics.active()) {
      neutronics.comm_.Bcast(densities_.data(), densities_.size(), MPI_DOUBLE);

      // For each OpenMC cell instance in a fluid cell, volume average the
      // densities and set
      // TODO:  Might be able to use xtensor masking to do some of this
      for (CellHandle cell = 0; cell < cell_to_elems_.size(); ++cell) {
        if (cell_fluid_mask_[cell] == 1) {
          double average_density = 0.0;
          double total_vol = 0.0;
          for (int e : cell_to_elems_.at(cell)) {
            average_density += densities_[e] * elem_volumes_[e];
            total_vol += elem_volumes_[e];
          }
          double density = average_density / total_vol;
          Ensures(density > 0.0);
          neutronics.set_density(cell, average_density / total_vol);
        }
      }
    }
  }
}

} // namespace enrico
