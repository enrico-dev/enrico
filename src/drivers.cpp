#include "drivers.h"
#include "mpi.h"
#include "nek_interface.h"
#include "openmc/capi.h"
#include "stream_geom.h"

#include <algorithm> // for max, fill, copy
#include <iterator> // for back_inserter
#include <unordered_set>

namespace stream {

// ============================================================================
// HeatFluids Driver
// ============================================================================

bool HeatFluidsDriver::active() const {
  return comm_.comm != MPI_COMM_NULL;
}

// ============================================================================
// Transport Driver
// ============================================================================

bool TransportDriver::active() const {
  return comm_.comm != MPI_COMM_NULL;
}

// ============================================================================
// OpenMC Driver
// ============================================================================

OpenmcDriver::OpenmcDriver(int argc, char *argv[], MPI_Comm comm)
    : TransportDriver(comm) {
  if (active()) {
    openmc_init(argc, argv, &comm);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void OpenmcDriver::init_step() { openmc_simulation_init(); }

void OpenmcDriver::solve_step() { openmc_run(); }

void OpenmcDriver::finalize_step() { openmc_simulation_finalize(); }

OpenmcDriver::~OpenmcDriver() {
  if (active())
    openmc_finalize();
  MPI_Barrier(MPI_COMM_WORLD);
}

// ============================================================================
// Nek5000 Driver
// ============================================================================

NekDriver::NekDriver(MPI_Comm comm) : HeatFluidsDriver(comm) {
  lelg_ = nek_get_lelg();
  lelt_ = nek_get_lelt();
  lx1_ = nek_get_lx1();

  if (active()) {
    MPI_Fint int_comm = MPI_Comm_c2f(comm_.comm);
    C2F_nek_init(static_cast<const int *>(&int_comm));

    nelgt_ = nek_get_nelgt();
    nelt_ = nek_get_nelt();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void NekDriver::init_step() {}

void NekDriver::solve_step() { C2F_nek_solve(); }

void NekDriver::finalize_step() {}

Position NekDriver::get_global_elem_centroid(int global_elem) const {
  Position centroid;
  int ierr = nek_get_global_elem_centroid(global_elem, &centroid);
  return centroid;
}

NekDriver::~NekDriver() {
  if (active())
    C2F_nek_end();
  MPI_Barrier(MPI_COMM_WORLD);
}

// ============================================================================
// OpenmcNekDriver
// ============================================================================

OpenmcNekDriver::OpenmcNekDriver(int argc, char **argv, MPI_Comm coupled_comm,
                                 MPI_Comm openmc_comm, MPI_Comm nek_comm, MPI_Comm intranode_comm) :
    comm_(coupled_comm),
    openmc_driver_(argc, argv, openmc_comm),
    nek_driver_(nek_comm),
    intranode_comm_(intranode_comm) {
  init_mappings();
  init_tallies();
};

void OpenmcNekDriver::init_mappings() {
  // TODO: This won't work if the Nek/OpenMC communicators are disjoint
  if (nek_driver_.active()) {
    // Create buffer to store material IDs corresponding to each Nek global
    // element. This is needed because calls to OpenMC API functions can only be
    // made from processes
    int32_t mat_idx[nek_driver_.nelgt_];

    if (openmc_driver_.active()) {
      std::unordered_set<int32_t> tracked;
      for (int global_elem = 1; global_elem <= nek_driver_.nelgt_; ++global_elem) {
        // Determine cell instance corresponding to global element
        Position elem_pos = nek_driver_.get_global_elem_centroid(global_elem);
        CellInstance c {elem_pos};
        if (tracked.find(c.material_index_) == tracked.end()) {
          openmc_driver_.cells_.push_back(c);
          tracked.insert(c.material_index_);
        }

        // Get corresponding material
        mat_to_elems_[c.material_index_].push_back(global_elem);

        // Set value for material ID in array
        mat_idx[global_elem - 1] = c.material_index_;
      }

      // Determine number of unique OpenMC materials
      n_materials_ = openmc_driver_.cells_.size();
    }

    // Broadcast array of material IDs to each Nek rank
    intranode_comm_.Bcast(mat_idx, nek_driver_.nelgt_, MPI_INT32_T);

    // Broadcast number of materials
    intranode_comm_.Bcast(&n_materials_, 1, MPI_INT32_T);

    // Set element -> material ID mapping on each Nek rank
    for (int global_elem = 1; global_elem <= nek_driver_.nelgt_; ++global_elem) {
      elem_to_mat_[global_elem] = mat_idx[global_elem - 1];
    }

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

void OpenmcNekDriver::init_tallies() {
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
    for (const auto &c : openmc_driver_.cells_) {
      mats.push_back(c.material_index_);
    }

    // Set bins for filter
    openmc_material_filter_set_bins(index_filter, mats.size(), mats.data());

    // Create tally and assign scores/filters
    openmc_extend_tallies(1, &openmc_driver_.index_tally_, nullptr);
    openmc_tally_allocate(openmc_driver_.index_tally_, "generic");
    openmc_tally_set_id(openmc_driver_.index_tally_, tally_id);
    char score_array[][20]{"kappa-fission"};
    const char *scores[]{score_array[0]}; // OpenMC expects a const char**, ugh
    openmc_tally_set_scores(openmc_driver_.index_tally_, 1, scores);
    openmc_tally_set_filters(openmc_driver_.index_tally_, 1, &index_filter);
  }
}

void OpenmcNekDriver::update_heat_source() {
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
      heat[i] = JOULE_PER_EV * results[3*i + 1] / m;

      // Sum up heat in each material
      total_heat += heat[i];
    }

    // TODO: Need to have total power in W specified by user
    double power = 1.0;

    // Normalize heat source in each material and collect in an array
    for (int i = 0; i < n_materials_; ++i) {
      // TODO: Need volumes from OpenMC
      double V = 1.0;

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

void OpenmcNekDriver::update_temperature() {
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
        displs[i] = displs[i-1] + recvcounts[i-1];
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
      std::copy(T, T+m, T_unordered);
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
