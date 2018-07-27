#include "drivers.h"
#include "mpi.h"
#include "nek_interface.h"
#include "openmc.h"
#include "stream_geom.h"

#include <algorithm> // for max
#include <unordered_set>

namespace stream {

// ============================================================================
// HeatFluids Driver
// ============================================================================

bool HeatFluidsDriver::active() const {
  return proc_info_.comm != MPI_COMM_NULL;
}

// ============================================================================
// Transport Driver
// ============================================================================

bool TransportDriver::active() const {
  return proc_info_.comm != MPI_COMM_NULL;
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

int32_t OpenmcDriver::get_mat_id(Position position) const {
  int32_t mat_id, instance;
  double xyz[3] = {position.x, position.y, position.z};
  openmc_find(xyz, 2, &mat_id, &instance);
  return mat_id;
}

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
    MPI_Fint int_comm = MPI_Comm_c2f(proc_info_.comm);
    C2F_nek_init(static_cast<const int *>(&int_comm));
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void NekDriver::init_step() {}

void NekDriver::solve_step() { C2F_nek_solve(); }

void NekDriver::finalize_step() {}

Position NekDriver::get_global_elem_centroid(int32_t global_elem) const {
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
    proc_info_(coupled_comm),
    openmc_driver_(argc, argv, openmc_comm),
    nek_driver_(nek_comm),
    intranode_(intranode_comm) {
  init_mappings();
  init_tallies();
};

void OpenmcNekDriver::init_mappings() {
  // TODO: This won't work if the Nek/OpenMC communicators are disjoint
  if (nek_driver_.active()) {
    // Create buffer to store material IDs corresponding to each Nek global
    // element. This is needed because calls to OpenMC API functions can only be
    // made from processes
    int32_t mat_ids[nek_driver_.lelg_];

    if (openmc_driver_.active()) {
      for (int global_elem = 1; global_elem <= nek_driver_.lelg_; ++global_elem) {
        Position elem_pos = nek_driver_.get_global_elem_centroid(global_elem);
        int32_t mat_id = openmc_driver_.get_mat_id(elem_pos);
        mats_to_elems_[mat_id].push_back(global_elem);

        // Set value for material ID in array
        mat_ids[global_elem - 1] = mat_id;
      }
    }

    // Broadcast array of material IDs to each Nek rank
    MPI_Bcast(mat_ids, nek_driver_.lelg_, MPI_INT32_T, 0, intranode_.comm);

    // Set element -> material ID mapping on each Nek rank
    for (int global_elem = 1; global_elem <= nek_driver_.lelg_; ++global_elem) {
      elems_to_mats_[global_elem] = mat_ids[global_elem - 1];
    }
  }
}

void OpenmcNekDriver::init_tallies() {
  if (openmc_driver_.active()) {
    // Determine maximum tally/filter ID used so far
    // TODO: Add functions in OpenMC API for this
    int32_t max_filter_id = 0;
    int32_t filter_id;
    for (int i = 1; i <= n_filters; ++i) {
      openmc_filter_get_id(i, &filter_id);
      max_filter_id = std::max(max_filter_id, filter_id);
    }

    int32_t max_tally_id = 0;
    int32_t tally_id;
    for (int i = 1; i <= n_tallies; ++i) {
      openmc_tally_get_id(i, &tally_id);
      max_tally_id = std::max(max_tally_id, tally_id);
    }

    int32_t index_filter;
    openmc_extend_filters(1, &index_filter, nullptr);
    openmc_filter_set_type(index_filter, "material");
    openmc_filter_set_id(index_filter, max_filter_id + 1);

    // Build set of material indices by looping over the OpenMC mat->Nek elem
    // mapping and using only the keys (material indices)
    std::unordered_set<int32_t> mat_set;
    for (const auto &pair : mats_to_elems_) {
      mat_set.insert(pair.first);
    }

    // Convert set to vector and then set bins for filter
    std::vector<int32_t> mats{mat_set.begin(), mat_set.end()};
    openmc_material_filter_set_bins(index_filter, mats.size(), mats.data());

    // Create tally and assign scores/filters
    openmc_extend_tallies(1, &openmc_driver_.index_tally_, nullptr);
    openmc_tally_set_type(openmc_driver_.index_tally_, "generic");
    openmc_tally_set_id(openmc_driver_.index_tally_, max_tally_id + 1);
    char score_array[][20]{"kappa-fission"};
    const char *scores[]{score_array[0]}; // OpenMC expects a const char**, ugh
    openmc_tally_set_scores(openmc_driver_.index_tally_, 1, scores);
    openmc_tally_set_filters(openmc_driver_.index_tally_, 1, &index_filter);
  }
}

void OpenmcNekDriver::update_heat_source() {

}

void OpenmcNekDriver::update_temperature() {

}

} // namespace stream
