#include "drivers.h"
#include "mpi.h"
#include "nek_interface.h"
#include "openmc.h"
#include "stream_geom.h"

#include <algorithm> // for max, fill
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
    int32_t mat_ids[nek_driver_.nelgt_];

    if (openmc_driver_.active()) {
      for (int global_elem = 1; global_elem <= nek_driver_.nelgt_; ++global_elem) {
        Position elem_pos = nek_driver_.get_global_elem_centroid(global_elem);
        int32_t mat_id = openmc_driver_.get_mat_id(elem_pos);
        mats_to_elems_[mat_id].push_back(global_elem);

        // Set value for material ID in array
        mat_ids[global_elem - 1] = mat_id;
      }

      // Determine number of unique OpenMC materials
      std::unordered_set<int32_t> mat_set;
      for (const auto& pair : mats_to_elems_) {
        mat_set.insert(pair.first);
      }
      n_materials_ = mat_set.size();
    }

    // Broadcast array of material IDs to each Nek rank
    intranode_comm_.Bcast(mat_ids, nek_driver_.nelgt_, MPI_INT32_T);

    // Broadcast number of materials
    intranode_comm_.Bcast(&n_materials_, 1, MPI_INT32_T);

    // Set element -> material ID mapping on each Nek rank
    for (int global_elem = 1; global_elem <= nek_driver_.nelgt_; ++global_elem) {
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

    int32_t& index_filter = openmc_driver_.index_filter_;
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
    // for each local element
    // get corresponding global element ID
    // if corresponding material
    // get corresponding material
    // get heat source for that material
    // Convert units from W/cm^3 to ???
    // set source for subsequent Nek run
  }
}

void OpenmcNekDriver::update_temperature() {
  if (nek_driver_.active()) {
    // Gather local temperatures into an array
    int n = nek_driver_.nelt_;
    double T_local[n];
    std::fill(T_local, T_local + n, 293.6);
    // TODO: Get temperature for each local element

    // Gather number of local elements into array on rank 0
    // TODO: Move this to initialization of Nek driver?
    int p = nek_driver_.comm_.rank == 0 ? nek_driver_.comm_.size : 0;
    int recvcounts[p];
    nek_driver_.comm_.Gather(&n, 1, MPI_INT, recvcounts, 1, MPI_INT);

    // Create array for all temperatures (size zero on non-root processes)
    int m = nek_driver_.comm_.rank == 0 ? nek_driver_.nelgt_ : 0;
    double T[m];

    // TODO: Need volumes of all global elements
    double vol[m];
    std::fill(vol, vol + m, 1.0);

    // collect temperatures from each local element onto root process
    if (nek_driver_.comm_.rank == 0) {
      // Create array of displacements for each process
      int displs[p];
      displs[0] = 0;
      for (int i = 1; i < p; ++i) {
        displs[i] = displs[i-1] + recvcounts[i-1];
      }

      // Gather temperatures onto root
      nek_driver_.comm_.Gatherv(
        T_local, n, MPI_DOUBLE,
        T, recvcounts, displs, MPI_DOUBLE
      );
    } else {
      // Send temperature to root via gather
      nek_driver_.comm_.Gatherv(
        T_local, n, MPI_DOUBLE,
        nullptr, nullptr, nullptr, MPI_DOUBLE
      );
    }
    // TODO: Temperatures need to be in order with respect to global element numbers

    if (openmc_driver_.active()) {
      // broadcast temperatures to all OpenMC processes
      openmc_driver_.comm_.Bcast(T, m, MPI_DOUBLE);
      openmc_driver_.comm_.Bcast(vol, m, MPI_DOUBLE);

      // For each OpenMC material, volume average temperatures and set
      for (const auto& pair : mats_to_elems_) {
        int32_t mat_id = pair.first;
        auto& global_elems = pair.second;

        // Get volume-average temperature for this material
        double average_temp = 0.0;
        for (int32_t elem : global_elems) {
          average_temp += T[elem - 1] * vol[elem - 1];
        }
        average_temp /= global_elems.size();

        // TODO: Set temperature
      }
    }
  }

}

} // namespace stream
