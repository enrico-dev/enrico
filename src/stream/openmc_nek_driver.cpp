#include "nek_interface.h"
#include "openmc/capi.h"
#include "openmc_nek_driver.h"
#include "stream_const.h"

namespace stream {

OpenmcNekDriver::OpenmcNekDriver(int argc, char** argv, MPI_Comm coupled_comm,
                                 MPI_Comm openmc_comm, MPI_Comm nek_comm, MPI_Comm intranode_comm) :
    comm_(coupled_comm),
    openmc_driver_(argc, argv, openmc_comm),
    nek_driver_(nek_comm),
    intranode_comm_(intranode_comm)
{
  init_mpi_datatypes();
  init_mappings();
  init_tallies();
};

void OpenmcNekDriver::init_mpi_datatypes() {
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

  int nlocals = nek_driver_.active() ? nek_driver_.nelt_ : 0;
  int nglobals = nek_driver_.active() ? nek_driver_.nelgt_ : 0;

  // Every Nek proc gets its local element centroids
  int lec_size = nek_driver_.active() ? nlocals : 0;
  std::vector<Position> local_element_centroids(lec_size);

  // Only the OpenMC procs get the global element centroids
  int gec_size = openmc_driver_.active() ? nglobals : 0;
  global_elem_centroids.reserve(gec_size);

  // Step 1: Get global element centroids on all OpenMC ranks
  if (nek_driver_.active()) {
    // Each Nek proc finds the centroids of its local elements
    for (int i = 0; i < nlocals; ++i) {
      local_element_centroids.at(i) = nek_driver_.get_local_elem_centroid(i);
    }
    // Gather all the local element centroids on the Nek5000/OpenMC root
    nek_driver_.comm_.Gatherv(local_element_centroids.data(), nlocals, position_mpi_datatype,
                              global_elem_centroids.data(), nek_driver_.local_counts_.data(),
                              nek_driver_.local_displs_.data(), position_mpi_datatype);
    // Broadcast global_element_centroids onto all the OpenMC procs
    if (openmc_driver_.active()) {
      openmc_driver_.comm_.Bcast(global_elem_centroids.data(), position_mpi_datatype);
    }
  }

  // Step 2: Set element->material and material->element mappings
  if (nek_driver_.active()) {
    // Create buffer to store material IDs corresponding to each Nek global
    // element. This is needed because calls to OpenMC API functions can only be
    // made from processes
    int32_t mat_idx[nglobals];

    if (openmc_driver_.active()) {
      std::unordered_set<int32_t> tracked;

      for (int i = 0; i < nglobals; ++i) {
        // Determine cell instance corresponding to global element
        Position elem_pos = global_elem_centroids[i];
        CellInstance c{elem_pos};
        if (tracked.find(c.material_index_) == tracked.end()) {
          openmc_driver_.cells_.push_back(c);
          tracked.insert(c.material_index_);
        }
        // Get corresponding material
        mat_to_elems_[c.material_index_].push_back(i);
        // Set value for material ID in array
        mat_idx[i - 1] = c.material_index_;
      }

      // Determine number of unique OpenMC materials
      n_materials_ = openmc_driver_.cells_.size();
    }

    // Broadcast array of material IDs to each Nek rank
    intranode_comm_.Bcast(mat_idx, nglobals, MPI_INT32_T);

    // Broadcast number of materials
    intranode_comm_.Bcast(&n_materials_, 1, MPI_INT32_T);

    // Set element -> material ID mapping on each Nek rank
    for (int i = 1; i <= nglobals; ++i) {
      elem_to_mat_[i] = mat_idx[i - 1];
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
    for (int local_elem = 1; local_elem <= nek_driver_.nelt_; ++local_elem) {
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
  if (nek_driver_.active()) {
    // Gather local temperatures into an array
    int n = nek_driver_.nelt_;
    double T_local[n];
    std::fill(T_local, T_local + n, 293.6);
    for (int local_elem = 1; local_elem <= n; ++local_elem) {
      // T_local[local_elem - 1] = nek_get_temperature(local_elem);
    }

    // Since local elements might be out of order with respect to global element
    // ordering, we need to know what global elements the local ones correspond
    // to
    int global_elems[n];
    for (int i = 0; i < n; ++i) {
      global_elems[i] = nek_get_global_elem(i + 1) - 1;
    }

    // Gather number of local elements into array on rank 0
    // TODO: Move this to initialization of Nek driver?
    int p = nek_driver_.comm_.rank == 0 ? nek_driver_.comm_.size : 0;
    int recvcounts[p];
    nek_driver_.comm_.Gather(&n, 1, MPI_INT, recvcounts, 1, MPI_INT);

    // Create array for all temperatures (size zero on non-root processes)
    int m = openmc_driver_.active() ? nek_driver_.nelgt_ : 0;
    double T[m];

    // Array for global element indices
    double all_global_elems[m];

    // TODO: Need volumes of all global elements
    double vol[m];
    std::fill(vol, vol + m, 1.0);
    for (int global_elem = 1; global_elem < m; ++global_elem) {
      // vol[global_elem - 1] = nek_get_volume(global_elem);
    }

    // collect temperatures from each local element onto root process
    if (nek_driver_.comm_.rank == 0) {
      // Create array of displacements for each process
      int displs[p];
      displs[0] = 0;
      for (int i = 1; i < p; ++i) {
        displs[i] = displs[i - 1] + recvcounts[i - 1];
      }

      // Gather temperatures and global element indices onto root
      nek_driver_.comm_.Gatherv(
          T_local, n, MPI_DOUBLE,
          T, recvcounts, displs, MPI_DOUBLE
      );
      nek_driver_.comm_.Gatherv(
          global_elems, n, MPI_INT,
          all_global_elems, recvcounts, displs, MPI_INT
      );
    } else {
      // Send temperature and global element indices to root via gather
      nek_driver_.comm_.Gatherv(
          T_local, n, MPI_DOUBLE,
          nullptr, nullptr, nullptr, MPI_DOUBLE
      );
      nek_driver_.comm_.Gatherv(
          global_elems, n, MPI_INT,
          nullptr, nullptr, nullptr, MPI_INT
      );
    }

    if (openmc_driver_.active()) {
      // broadcast data from Nek root to all OpenMC processes
      openmc_driver_.comm_.Bcast(T, m, MPI_DOUBLE);
      openmc_driver_.comm_.Bcast(all_global_elems, m, MPI_INT);
      openmc_driver_.comm_.Bcast(vol, m, MPI_DOUBLE);

      // At this point, the values in T are not in the correct ordering with
      // respect to the global element array. Here, we make a copy of T and
      // reorder them according to the position they should be in
      // (all_global_elems)
      double T_unordered[m];
      std::copy(T, T + m, T_unordered);
      for (int i = 0; i < m; ++i) {
        int global_index = all_global_elems[i];
        T[global_index] = T_unordered[i];
      }

      // For each OpenMC material, volume average temperatures and set
      for (const auto& c : openmc_driver_.cells_) {
        // Get corresponding global elements
        const auto& global_elems = mat_to_elems_.at(c.material_index_);

        // Get volume-average temperature for this material
        double average_temp = 0.0;
        double total_vol = 0.0;
        for (int elem : global_elems) {
          average_temp += T[elem - 1] * vol[elem - 1];
          total_vol += vol[elem - 1];
        }
        average_temp /= total_vol;

        // Set temperature for cell instance
        c.set_temperature(average_temp);
      }
    }
  }

}

} // namespace stream
