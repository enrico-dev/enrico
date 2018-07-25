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
  lelg = nek_get_lelg();
  lelt = nek_get_lelt();
  lx1 = nek_get_lx1();

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

OpenmcNekDriver::OpenmcNekDriver(int argc, char **argv, MPI_Comm coupled_comm, MPI_Comm openmc_comm, MPI_Comm nek_comm) :
    proc_info_(coupled_comm),
    openmc_driver_(argc, argv, openmc_comm),
    nek_driver_(nek_comm) {
  init_mats_to_elems();
  init_tallies();
};

void OpenmcNekDriver::init_mats_to_elems() {
  if (openmc_driver_.active()) {
    for (int global_elem = 1; global_elem <= nek_driver_.lelg; ++global_elem) {
      Position elem_pos = nek_driver_.get_global_elem_centroid(global_elem);
      int32_t mat_id = openmc_driver_.get_mat_id(elem_pos);
      mats_to_elems_[mat_id].push_back(global_elem);
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
    openmc_extend_tallies(1, &openmc_driver_.index_tally, nullptr);
    openmc_tally_set_type(openmc_driver_.index_tally, "generic");
    openmc_tally_set_id(openmc_driver_.index_tally, max_tally_id + 1);
    char score_array[][20]{"kappa-fission"};
    const char *scores[]{score_array[0]}; // OpenMC expects a const char**, ugh
    openmc_tally_set_scores(openmc_driver_.index_tally, 1, scores);
    openmc_tally_set_filters(openmc_driver_.index_tally, 1, &index_filter);
  }
}

} // namespace stream
