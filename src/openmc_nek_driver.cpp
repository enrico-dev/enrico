#include "enrico/openmc_nek_driver.h"

#include "enrico/const.h"
#include "enrico/error.h"
#include "enrico/nek_interface.h"
#include "enrico/message_passing.h"

#include "heat_xfer_backend.h"
#include "openmc/capi.h"
#include "pugixml.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xnorm.hpp"
#include <gsl/gsl>

#include <string>

namespace enrico {

OpenmcNekDriver::OpenmcNekDriver(MPI_Comm comm, pugi::xml_node node) :
  CoupledDriver{comm, node}
{
  // Get parameters from enrico.xml
  pugi::xml_node nek_node = node.child("nek5000");
  pressure_ = node.child("pressure").text().as_double();
  openmc_procs_per_node_ = node.child("openmc_procs_per_node").text().as_int();

  // Postcondition checks on user inputs
  Expects(pressure_ > 0.0);
  Expects(openmc_procs_per_node_ > 0);

  // Create communicator for OpenMC with 1 process per node
  MPI_Comm openmc_comm;
  MPI_Comm intranode_comm;
  enrico::get_node_comms(comm_.comm, openmc_procs_per_node_, &openmc_comm, &intranode_comm);

  // Set intranode communicator
  intranode_comm_ = Comm(intranode_comm);

  // Instantiate OpenMC and Nek drivers
  openmc_driver_ = std::make_unique<OpenmcDriver>(openmc_comm);
  nek_driver_ = std::make_unique<NekDriver>(comm, nek_node);

  // Determine number of local/global elements for each rank
  n_local_elem_ = nek_driver_->active() ? nek_driver_->nelt_ : 0;
  n_global_elem_ = nek_driver_->active() ? nek_driver_->nelgt_ : 0;

  init_mpi_datatypes();
  init_mappings();
  init_tallies();
  init_volumes();
  init_temperatures();
  init_densities();
};

OpenmcNekDriver::~OpenmcNekDriver()
{
  free_mpi_datatypes();
}

Driver& OpenmcNekDriver::getNeutronicsDriver() const
{
  return *openmc_driver_;
}

Driver& OpenmcNekDriver::getHeatDriver() const
{
  return *nek_driver_;
}

void OpenmcNekDriver::init_mpi_datatypes()
{
  // Currently, this sets up only position_mpi_datatype
  Position p;
  int blockcounts[3] = {1, 1, 1};
  MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint displs[3];

  // Get displacements of struct members
  MPI_Get_address(&p.x, &displs[0]);
  MPI_Get_address(&p.y, &displs[1]);
  MPI_Get_address(&p.z, &displs[2]);

  // Make the displacements relative
  displs[2] -= displs[0];
  displs[1] -= displs[0];
  displs[0] = 0;

  // Make datatype
  MPI_Type_create_struct(3, blockcounts, displs, types, &position_mpi_datatype);
  MPI_Type_commit(&position_mpi_datatype);
}

void OpenmcNekDriver::init_mappings()
{
  // TODO: This won't work if the Nek/OpenMC communicators are disjoint

  if (nek_driver_->active()) {
    // Only the OpenMC procs get the global element centroids/fluid-identities
    if (openmc_driver_->active()) {
      elem_centroids_.resize(n_global_elem_);
      elem_is_in_fluid_.resize(n_global_elem_);
    }
    // Step 1: Get global element centroids/fluid-identities on all OpenMC ranks
    // Each Nek proc finds the centroids/fluid-identities of its local elements
    Position local_element_centroids[n_local_elem_];
    int local_element_is_in_fluid[n_local_elem_];
    for (int i = 0; i < n_local_elem_; ++i) {
      local_element_centroids[i] = nek_driver_->get_local_elem_centroid(i+1);
      local_element_is_in_fluid[i] = nek_driver_->local_elem_is_in_fluid(i+1);
    }
    // Gather all the local element centroids/fluid-identities on the Nek5000/OpenMC root
    nek_driver_->comm_.Gatherv(local_element_centroids, n_local_elem_, position_mpi_datatype,
                              elem_centroids_.data(), nek_driver_->local_counts_.data(),
                              nek_driver_->local_displs_.data(), position_mpi_datatype);
    nek_driver_->comm_.Gatherv(local_element_is_in_fluid, n_local_elem_, MPI_INT,
                               elem_is_in_fluid_.data(), nek_driver_->local_counts_.data(),
                               nek_driver_->local_displs_.data(), MPI_INT);
    // Broadcast global_element_centroids/fluid-identities onto all the OpenMC procs
    if (openmc_driver_->active()) {
      openmc_driver_->comm_.Bcast(elem_centroids_.data(), n_global_elem_,
                                 position_mpi_datatype);
      openmc_driver_->comm_.Bcast(elem_is_in_fluid_.data(), n_global_elem_, MPI_INT);
    }

    // Step 2: Set element->material and material->element mappings
    // Create buffer to store material IDs corresponding to each Nek global
    // element. This is needed because calls to OpenMC API functions can only be
    // made from processes
    int32_t mat_idx[n_global_elem_];

    if (openmc_driver_->active()) {
      std::unordered_set<int32_t> tracked;

      for (int i = 0; i < n_global_elem_; ++i) {
        // Determine cell instance corresponding to global element
        Position elem_pos = elem_centroids_[i];
        CellInstance c {elem_pos};
        if (tracked.find(c.material_index_) == tracked.end()) {
          openmc_driver_->cells_.push_back(c);
          tracked.insert(c.material_index_);
        }
        // Get corresponding material
        mat_to_elems_[c.material_index_].push_back(i);
        // Set value for material ID in array
        mat_idx[i] = c.material_index_;
      }

      // Determine number of unique OpenMC materials
      n_materials_ = openmc_driver_->cells_.size();
    }

    // Broadcast array of material IDs to each Nek rank
    intranode_comm_.Bcast(mat_idx, n_global_elem_, MPI_INT32_T);

    // Broadcast number of materials
    intranode_comm_.Bcast(&n_materials_, 1, MPI_INT32_T);

    // Set element -> material ID mapping on each Nek rank
    for (int i = 0; i < n_global_elem_; ++i) {
      elem_to_mat_[i] = mat_idx[i];
    }

    // Step 3: Map material indices to positions in heat array
    // Each Nek rank needs to know where to find heat source values in an array.
    // This requires knowing a mapping of material indices to positions in the
    // heat array (of size n_materials_). For this mapping, we use a direct
    // address table and begin by filling it with zeros.
    heat_index_.reserve(n_materials_);
    std::fill_n(std::back_inserter(heat_index_), n_materials_, 0);

    // Now, set the mapping values on the OpenMC processes...
    if (openmc_driver_->active()) {
      for (int i = 0; i < n_materials_; ++i) {
        int32_t mat_index = openmc_driver_->cells_[i].material_index_;
        heat_index_.at(mat_index) = i;
      }
    }

    // ...and broadcast to the other Nek processes
    intranode_comm_.Bcast(heat_index_.data(), n_materials_, MPI_INT);
  }
}

void OpenmcNekDriver::init_tallies()
{
  if (openmc_driver_->active()) {
    // Build vector of material indices
    std::vector<int32_t> mats;
    for (const auto& c : openmc_driver_->cells_) {
      mats.push_back(c.material_index_);
    }
    openmc_driver_->create_tallies(mats);
  }
}

void OpenmcNekDriver::init_temperatures()
{
  // TODO: This won't work if the Nek/OpenMC communicators are disjoint
  // Only the OpenMC procs get the global temperatures
  if (nek_driver_->active() and openmc_driver_->active()) {
    elem_temperatures_.resize({n_global_elem_});
    elem_temperatures_prev_.resize({n_global_elem_});

    std::fill(elem_temperatures_.begin(), elem_temperatures_.end(), 293.6);
    std::fill(elem_temperatures_prev_.begin(), elem_temperatures_prev_.end(), 293.6);
  }
}

void OpenmcNekDriver::init_volumes()
{
  // TODO: This won't work if the Nek/OpenMC communicators are disjoint

  // Only the OpenMC procs get the global element volumes (gev)
  if (nek_driver_->active()) {
    if (openmc_driver_->active()) {
      elem_volumes_.resize({n_global_elem_});
    }
    // Every Nek proc gets its local element volumes (lev)
    double local_elem_volumes[n_local_elem_];
    for (int i = 0; i < n_local_elem_; ++i) {
      local_elem_volumes[i] = nek_driver_->get_local_elem_volume(i + 1);
    }
    // Gather all the local element volumes on the Nek5000/OpenMC root
    nek_driver_->comm_.Gatherv(local_elem_volumes, n_local_elem_, MPI_DOUBLE,
                              elem_volumes_.data(), nek_driver_->local_counts_.data(),
                              nek_driver_->local_displs_.data(), MPI_DOUBLE);
    // Broadcast global_element_volumes onto all the OpenMC procs
    if (openmc_driver_->active()) {
      openmc_driver_->comm_.Bcast(elem_volumes_.data(), n_global_elem_, MPI_DOUBLE);
    }
  }
}

void OpenmcNekDriver::init_densities()
{
  // TODO: This won't work if the Nek/OpenMC communicators are disjoint

  // Only the OpenMC procs get space for the global densities
  // Note that the values of the densities are not initialized!
  if (nek_driver_->active() and openmc_driver_->active()) {
    elem_densities_.resize({n_global_elem_});
  }
}

void OpenmcNekDriver::update_heat_source()
{
  // Create array to store volumetric heat deposition in each material
  xt::xtensor<double, 1> heat = xt::empty<double>({n_materials_});

  if (openmc_driver_->active()) {
    // Get heat source normalized by user-specified power
    heat = openmc_driver_->heat_source(power_);
  }

  // OpenMC has heat source on each of its ranks. We need to make heat
  // source available on each Nek rank.
  intranode_comm_.Bcast(heat.data(), n_materials_, MPI_DOUBLE);

  if (nek_driver_->active()) {
    // Determine displacement for this rank
    int displacement = nek_driver_->local_displs_[nek_driver_->comm_.rank];

    // Loop over local elements to set heat source
    for (int local_elem = 1; local_elem <= n_local_elem_; ++local_elem) {
      // get corresponding global element
      int global_index = local_elem + displacement - 1;

      // get corresponding material
      int32_t mat_index = elem_to_mat_.at(global_index);
      int i = heat_index_.at(mat_index);

      err_chk(nek_set_heat_source(local_elem, heat[i]),
          "Error setting heat source for local element " + std::to_string(i));
    }
  }
}

void OpenmcNekDriver::update_temperature()
{
  if (nek_driver_->active()) {
    if (openmc_driver_->active()) {
      std::copy(elem_temperatures_.begin(), elem_temperatures_.end(),
          elem_temperatures_prev_.begin());
    }
    // Each Nek proc finds the temperatures of its local elements
    double local_elem_temperatures[n_local_elem_];
    for (int i = 0; i < n_local_elem_; ++i) {
      local_elem_temperatures[i] = nek_driver_->get_local_elem_temperature(i + 1);
    }
    // Gather all the local element temperatures on the Nek5000/OpenMC root
    nek_driver_->comm_.Gatherv(local_elem_temperatures, n_local_elem_, MPI_DOUBLE,
                              elem_temperatures_.data(), nek_driver_->local_counts_.data(),
                              nek_driver_->local_displs_.data(), MPI_DOUBLE);

    if (openmc_driver_->active()) {
      // Broadcast global_element_temperatures onto all the OpenMC procs
      openmc_driver_->comm_.Bcast(elem_temperatures_.data(), n_global_elem_, MPI_DOUBLE);

      // For each OpenMC material, volume average temperatures and set
      for (const auto& c : openmc_driver_->cells_) {

        // Get corresponding global elements
        const auto& global_elems = mat_to_elems_.at(c.material_index_);

        // Get volume-average temperature for this material
        double average_temp = 0.0;
        double total_vol = 0.0;
        for (int elem : global_elems) {
          double T = elem_temperatures_[elem];
          double V = elem_volumes_[elem];
          average_temp += T*V;
          total_vol += V;
        }

        // Set temperature for cell instance
        average_temp /= total_vol;
        c.set_temperature(average_temp);
      }
    }
  }
}

void OpenmcNekDriver::update_density()
{
  if (nek_driver_->active())
  {
    if (openmc_driver_->active())
    {
      for (const auto& c: openmc_driver_->cells_)
      {
        // Get corresponding global elements
        const auto& global_elems = mat_to_elems_.at(c.material_index_);

        bool any_in_fluid = false;
        for (int elem : global_elems)
        {
          if (elem_is_in_fluid_[elem] == 1)
          {
            any_in_fluid = true;
            break;
          }
        }

        double average_density = 0.0;
        double total_vol = 0.0;
        for (int elem : global_elems)
        {
          double T = elem_temperatures_[elem];
          double V = elem_volumes_[elem];

          if (any_in_fluid)
          {
            // nu1 returns specific volume in [m^3/kg]
            double density = 1.0e-3/nu1(pressure_, T);
            average_density += density*V;
            total_vol += V;
          }
        }

        // Set density for cell instance
        if (any_in_fluid)
          c.set_density(average_density / total_vol);
      }
    }
  }
}

bool OpenmcNekDriver::is_converged()
{
  bool converged;
  // WARNING: Assumes that OpenmcNekDriver rank 0 is in openmc_driver_->comm
  if (comm_.rank == 0) {
    auto n_expr = xt::norm_linf(elem_temperatures_ - elem_temperatures_prev_);
    auto n_scalar = n_expr();
    converged = (n_scalar < epsilon_);

    std::string msg = "temperature norm_linf: " + std::to_string(n_scalar);
    comm_.message(msg);
  }
  err_chk(comm_.Bcast(&converged, 1, MPI_CXX_BOOL));
  return converged;
}

void OpenmcNekDriver::free_mpi_datatypes()
{
  MPI_Type_free(&position_mpi_datatype);
}

} // namespace enrico
