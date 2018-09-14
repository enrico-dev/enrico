#include "nek_interface.h"
#include "openmc/capi.h"
#include "openmc_nek_driver.h"
#include "stream_const.h"
#include <iostream>

namespace stream {

OpenmcNekDriver::OpenmcNekDriver(int argc, char** argv, MPI_Comm coupled_comm,
                                 MPI_Comm openmc_comm, MPI_Comm nek_comm, MPI_Comm intranode_comm) :
    comm_(coupled_comm),
    openmc_driver_(argc, argv, openmc_comm),
    nek_driver_(nek_comm),
    intranode_comm_(intranode_comm)
{
  local_vector_size_ = nek_driver_.active() ? nek_driver_.nelt_ : 0;
  global_vector_size_ = openmc_driver_.active() ? nek_driver_.nelgt_ : 0;

  init_mpi_datatypes();
  init_mappings();
  init_tallies();
  init_volumes();
  init_temperatures();
};

OpenmcNekDriver::~OpenmcNekDriver()
{
  free_mpi_datatypes();
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
  displs[2] -= displs[1];
  displs[1] -= displs[0];
  displs[0] = 0;

  // Make datatype
  MPI_Type_create_struct(3, blockcounts, displs, types, &position_mpi_datatype);
  MPI_Type_commit(&position_mpi_datatype);
}

void OpenmcNekDriver::init_mappings()
{
  // TODO: This won't work if the Nek/OpenMC communicators are disjoint

  // Every Nek proc gets its local element centroids
  Position local_element_centroids[local_vector_size_];

  // Only the OpenMC procs get the global element centroids
  global_elem_centroids_.resize(global_vector_size_);

  // Step 1: Get global element centroids on all OpenMC ranks
  if (nek_driver_.active()) {
    // Each Nek proc finds the centroids of its local elements
    for (int i = 0; i < local_vector_size_; ++i) {
      local_element_centroids[i] = nek_driver_.get_local_elem_centroid(i);
    }
    // Gather all the local element centroids on the Nek5000/OpenMC root
    nek_driver_.comm_.Gatherv(local_element_centroids, local_vector_size_, position_mpi_datatype,
                              global_elem_centroids_.data(), nek_driver_.local_counts_.data(),
                              nek_driver_.local_displs_.data(), position_mpi_datatype);
    // Broadcast global_element_centroids onto all the OpenMC procs
    if (openmc_driver_.active()) {
      openmc_driver_.comm_.Bcast(global_elem_centroids_.data(), global_vector_size_,
                                 position_mpi_datatype);
    }
  }

  // Step 2: Set element->material and material->element mappings
  if (nek_driver_.active()) {
    // Create buffer to store material IDs corresponding to each Nek global
    // element. This is needed because calls to OpenMC API functions can only be
    // made from processes
    int mat_idx_size = nek_driver_.nelgt_;
    int32_t mat_idx[mat_idx_size];

    if (openmc_driver_.active()) {
      std::unordered_set<int32_t> tracked;

      for (int i = 0; i < mat_idx_size; ++i) {
        // Determine cell instance corresponding to global element
        Position elem_pos = global_elem_centroids_[i];
        CellInstance c{elem_pos};
        if (tracked.find(c.material_index_) == tracked.end()) {
          openmc_driver_.cells_.push_back(c);
          tracked.insert(c.material_index_);
        }
        // Get corresponding material
        mat_to_elems_[c.material_index_].push_back(i);
        // Set value for material ID in array
        mat_idx[i] = c.material_index_;
      }

      // Determine number of unique OpenMC materials
      n_materials_ = openmc_driver_.cells_.size();
    }

    // Broadcast array of material IDs to each Nek rank
    intranode_comm_.Bcast(mat_idx, mat_idx_size, MPI_INT32_T);

    // Broadcast number of materials
    intranode_comm_.Bcast(&n_materials_, 1, MPI_INT32_T);

    // Set element -> material ID mapping on each Nek rank
    for (int i = 1; i <= mat_idx_size; ++i) {
      elem_to_mat_[i] = mat_idx[i];
    }
  }

  // Step 3: Map material indices to positions in heat array
  if (nek_driver_.active()) {
    // Each Nek rank needs to know where to find heat source values in an array.
    // This requires knowing a mapping of material indices to positions in the
    // heat array (of size n_materials_). For this mapping, we use a direct
    // address table and begin by filling it with zeros.
    heat_index_.reserve(n_materials_);
    std::fill_n(std::back_inserter(heat_index_), n_materials_, 0);

    // Now, set the mapping values on the OpenMC processes...
    if (openmc_driver_.active()) {
      for (int i = 0; i < n_materials_; ++i) {
        int32_t mat_index = openmc_driver_.cells_[i].material_index_ - 1;
        heat_index_.at(mat_index) = i;
      }
    }

    // ...and broadcast to the other Nek processes
    intranode_comm_.Bcast(heat_index_.data(), n_materials_, MPI_INT);
  }
}

void OpenmcNekDriver::init_tallies()
{
  if (openmc_driver_.active()) {
    // Determine maximum tally/filter ID used so far
    int32_t filter_id, tally_id;
    openmc_get_filter_next_id(&filter_id);
    openmc_get_tally_next_id(&tally_id);

    int32_t& index_filter = openmc_driver_.index_filter_;
    openmc_extend_filters(1, &index_filter, nullptr);
    openmc_filter_set_type(index_filter, "material");
    openmc_filter_set_id(index_filter, filter_id);

    // Build vector of material indices
    std::vector<int32_t> mats;
    for (const auto& c : openmc_driver_.cells_) {
      mats.push_back(c.material_index_);
    }

    // Set bins for filter
    openmc_material_filter_set_bins(index_filter, mats.size(), mats.data());

    // Create tally and assign scores/filters
    openmc_extend_tallies(1, &openmc_driver_.index_tally_, nullptr);
    openmc_tally_allocate(openmc_driver_.index_tally_, "generic");
    openmc_tally_set_id(openmc_driver_.index_tally_, tally_id);
    char score_array[][20]{"kappa-fission"};
    const char* scores[]{score_array[0]}; // OpenMC expects a const char**, ugh
    openmc_tally_set_scores(openmc_driver_.index_tally_, 1, scores);
    openmc_tally_set_filters(openmc_driver_.index_tally_, 1, &index_filter);
  }
}

void OpenmcNekDriver::init_temperatures()
{
  // TODO: This won't work if the Nek/OpenMC communicators are disjoint
  // Only the OpenMC procs get the global temperatures
  global_elem_temperatures_.resize(global_vector_size_);
  std::fill(global_elem_temperatures_.begin(), global_elem_temperatures_.end(), 293.6);
}

void OpenmcNekDriver::init_volumes()
{
  // TODO: This won't work if the Nek/OpenMC communicators are disjoint

  // Every Nek proc gets its local element volumes (lev)
  double local_elem_volumes[local_vector_size_];

  // Only the OpenMC procs get the global element volumes (gev)
  global_elem_volumes_.resize(global_vector_size_);

  if (nek_driver_.active()) {
    // Each Nek proc finds the volumes of its local elements
    for (int i = 0; i < local_vector_size_; ++i) {
      local_elem_volumes[i] = nek_driver_.get_local_elem_volume(i + 1);
    }
    // Gather all the local element volumes on the Nek5000/OpenMC root
    nek_driver_.comm_.Gatherv(local_elem_volumes, local_vector_size_, MPI_DOUBLE,
                              global_elem_volumes_.data(), nek_driver_.local_counts_.data(),
                              nek_driver_.local_displs_.data(), MPI_DOUBLE);
    // Broadcast global_element_volumes onto all the OpenMC procs
    if (openmc_driver_.active()) {
      openmc_driver_.comm_.Bcast(global_elem_volumes_.data(), global_vector_size_, MPI_DOUBLE);
    }
  }
}

void OpenmcNekDriver::update_heat_source()
{
  // Create array to store volumetric heat deposition in each material
  double heat[n_materials_];

  if (openmc_driver_.active()) {
    // Get material bins
    int32_t* mats;
    int32_t n_mats;
    openmc_material_filter_get_bins(openmc_driver_.index_filter_, &mats, &n_mats);

    // Get tally results and number of realizations
    double* results;
    int shape[3];
    openmc_tally_results(openmc_driver_.index_tally_, &results, shape);
    int32_t m;
    openmc_tally_get_n_realizations(openmc_driver_.index_tally_, &m);

    // Determine energy production in each material
    double total_heat = 0.0;
    for (int i = 0; i < n_materials_; ++i) {
      // Get mean value for tally result and convert units from eV to J
      // TODO: Get rid of flattened array index by using xtensor?
      heat[i] = JOULE_PER_EV * results[3 * i + 1] / m;

      // Sum up heat in each material
      total_heat += heat[i];
    }

    // TODO: Need to have total power in W specified by user
    double power = 1.0;

    // Normalize heat source in each material and collect in an array
    for (int i = 0; i < n_materials_; ++i) {
      // Get volume
      double V = openmc_driver_.cells_.at(i).volume_;

      // Convert heat from J/src to W/cm^3. Dividing by total_heat gives the
      // fraction of heat deposited in each material. Multiplying by power
      // givens an absolute value in W
      double normalization = power / (total_heat * V);
      heat[i] *= normalization;
    }
  }

  // OpenMC has heat source on each of its ranks. We need to make heat
  // source available on each Nek rank.
  intranode_comm_.Bcast(heat, n_materials_, MPI_DOUBLE);

  if (nek_driver_.active()) {
    // Loop over local elements to set heat source
    for (int local_elem = 1; local_elem <= local_vector_size_; ++local_elem) {
      // get corresponding global element ID
      int global_elem = nek_get_global_elem(local_elem);

      // get corresponding material
      int32_t mat_index = elem_to_mat_.at(global_elem);
      int i = get_heat_index(mat_index);

      // TODO: set source for subsequent Nek run
      // nek_set_heat_source(local_elem, heat[i]);
    }
  }
}

void OpenmcNekDriver::update_temperature()
{
  // Every Nek proc gets its local element temperatures
  double local_elem_temperatures[local_vector_size_];

  if (nek_driver_.active()) {
    // Each Nek proc finds the temperatures of its local elements
    for (int i = 0; i < local_vector_size_; ++i) {
      local_elem_temperatures[i] = nek_driver_.get_local_elem_temperature(i + 1);
    }
    // Gather all the local element temperatures on the Nek5000/OpenMC root
    nek_driver_.comm_.Gatherv(local_elem_temperatures, local_vector_size_, MPI_DOUBLE,
                              global_elem_temperatures_.data(), nek_driver_.local_counts_.data(),
                              nek_driver_.local_displs_.data(), MPI_DOUBLE);

    if (openmc_driver_.active()) {
      // Broadcast global_element_temperatures onto all the OpenMC procs
      openmc_driver_.comm_.Bcast(global_elem_temperatures_.data(), global_vector_size_, MPI_DOUBLE);

      // For each OpenMC material, volume average temperatures and set
      for (const auto& c : openmc_driver_.cells_) {

        // Get corresponding global elements
        const auto& global_elems = mat_to_elems_.at(c.material_index_);

        // Get volume-average temperature for this material
        double average_temp = 0.0;
        double total_vol = 0.0;
        for (int elem : global_elems) {
          average_temp += global_elem_temperatures_[elem] * global_elem_volumes_[elem];
          total_vol += global_elem_volumes_[elem];
        }
        average_temp /= total_vol;

        // Set temperature for cell instance
        c.set_temperature(average_temp);
      }
    }
  }
}

void OpenmcNekDriver::free_mpi_datatypes()
{
  MPI_Type_free(&position_mpi_datatype);
}

} // namespace stream
