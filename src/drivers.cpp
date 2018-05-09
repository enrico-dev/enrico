#include "drivers.h"
#include "mpi.h"
#include "openmc.h"
#include "nek_interface.h"

#include <unordered_set>

// ============================================================================
// OpenMC Driver
// ============================================================================

OpenmcDriver::OpenmcDriver(MPI_Comm comm) : NeutronDriver(comm) {
  if (procInfo.comm != MPI_COMM_NULL) {
    // ROR: 2018-03-22: MPI_Comm_c2f is a macro (in MPICH, at least),
    // so we can't pass something like:
    //     openmc_init(&MPI_Comm_c2f(comm));
    // Hence, the dummy variable.
    MPI_Fint intComm = MPI_Comm_c2f(procInfo.comm);
    openmc_init(static_cast<const int *>(&intComm));

    int32_t index_filters[1];
    openmc_extend_filters(1, &index_filters[0], nullptr);
    openmc_filter_set_type(index_filters[0], "material");

    // Build set of material indices by looping over the OpenMC mat->Nek elem
    // mapping and using only the keys (material indices)
    std::unordered_set<int32_t> mat_set;
    for (const auto &pair : mats_to_nek_elems) {
      mat_set.insert(pair.first);
    }

    // Convert set to vector and then set bins for filter
    std::vector<int32_t> mats{mat_set.begin(), mat_set.end()};
    openmc_material_filter_set_bins(index_filters[0], mats.size(), mats.data());

    // Create tally and assign scores/filters
    int32_t index_tally;
    openmc_extend_tallies(1, &index_tally, nullptr);
    openmc_tally_set_type(index_tally, "generic");
    char score_array[][20]{"kappa-fission"};
    const char *scores[]{score_array[0]}; // OpenMC expects a const char**, ugh
    openmc_tally_set_scores(index_tally, 1, scores);
    openmc_tally_set_filters(index_tally, 1, index_filters);

    openmc_simulation_init();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void OpenmcDriver::initStep() {
}

void OpenmcDriver::solveStep() {
  if (procInfo.comm != MPI_COMM_NULL)
    openmc_run();
  MPI_Barrier(MPI_COMM_WORLD);
}

void OpenmcDriver::finalizeStep() {
  if (procInfo.comm != MPI_COMM_NULL)
    openmc_simulation_finalize();
  MPI_Barrier(MPI_COMM_WORLD);
}

OpenmcDriver::~OpenmcDriver() {
  if (procInfo.comm != MPI_COMM_NULL)
    openmc_finalize();
  MPI_Barrier(MPI_COMM_WORLD);
}

// ============================================================================
// Nek5000 Driver
// ============================================================================

NekDriver::NekDriver(MPI_Comm comm) : ThDriver(comm) {
  // ROR: 2018-03-22: MPI_Comm_c2f is a macro (in MPICH, at least),
  // so we can't pass something like:
  //     openmc_init(&MPI_Comm_c2f(comm));
  // Hence, the dummy variable.
  if (procInfo.comm != MPI_COMM_NULL) {
    MPI_Fint intComm = MPI_Comm_c2f(procInfo.comm);
    C2F_nek_init(static_cast<const int *>(&intComm));
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void NekDriver::initStep() {
}

void NekDriver::solveStep() {
  if (procInfo.comm != MPI_COMM_NULL)
    C2F_nek_solve();
  MPI_Barrier(MPI_COMM_WORLD);
}

void NekDriver::finalizeStep() {
}

NekDriver::~NekDriver() {
  if (procInfo.comm != MPI_COMM_NULL)
    C2F_nek_end();
  MPI_Barrier(MPI_COMM_WORLD);
}
