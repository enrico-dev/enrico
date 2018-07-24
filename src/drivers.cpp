#include "drivers.h"
#include "mpi.h"
#include "nek_interface.h"
#include "openmc.h"
#include "stream_geom.h"

#include <algorithm> // for max
#include <unordered_set>

// ============================================================================
// OpenMC Driver
// ============================================================================

OpenmcDriver::OpenmcDriver(int argc, char *argv[], MPI_Comm comm)
    : NeutronDriver(comm) {
  if (procInfo.comm != MPI_COMM_NULL) {
    openmc_init(argc, argv, &comm);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void OpenmcDriver::initStep() { openmc_simulation_init(); }

void OpenmcDriver::solveStep() { openmc_run(); }

void OpenmcDriver::finalizeStep() { openmc_simulation_finalize(); }

int32_t OpenmcDriver::getMatId(Position position) const {
  int32_t matId, instance;
  double xyz[3] = {position.x, position.y, position.z};
  openmc_find(xyz, 2, &matId, &instance);
  return matId;
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
  lelg = nek_get_lelg();
  lelt = nek_get_lelt();
  lx1 = nek_get_lx1();

  if (procInfo.comm != MPI_COMM_NULL) {
    MPI_Fint intComm = MPI_Comm_c2f(procInfo.comm);
    C2F_nek_init(static_cast<const int *>(&intComm));
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void NekDriver::initStep() {}

void NekDriver::solveStep() { C2F_nek_solve(); }

void NekDriver::finalizeStep() {}

Position NekDriver::getGlobalElemCentroid(const int32_t globalElem) {
  Position centroid;
  int ierr = nek_get_global_elem_centroid(globalElem, &centroid);
  return centroid;
}

NekDriver::~NekDriver() {
  if (procInfo.comm != MPI_COMM_NULL)
    C2F_nek_end();
  MPI_Barrier(MPI_COMM_WORLD);
}

// ============================================================================
// OpenmcNekDriver
// ============================================================================

OpenmcNekDriver::OpenmcNekDriver(int argc, char **argv, MPI_Comm coupledComm, MPI_Comm openmcComm, MPI_Comm nekComm) :
    procInfo(coupledComm),
    openmcDriver(argc, argv, openmcComm),
    nekDriver(nekComm) {
  initMatsToElems();
  initTallies();
};

void OpenmcNekDriver::initMatsToElems() {
  if (openmcDriver.procInfo.comm != MPI_COMM_NULL) {
    for (int globalElem = 1; globalElem <= nekDriver.lelg; ++globalElem) {
      Position elemPos = nekDriver.getGlobalElemCentroid(globalElem);
      int32_t matId = openmcDriver.getMatId(elemPos);
      matsToElems[matId].push_back(globalElem);
    }
  }
}

void OpenmcNekDriver::initTallies() {
  if (openmcDriver.procInfo.comm != MPI_COMM_NULL) {
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
    for (const auto &pair : matsToElems) {
      mat_set.insert(pair.first);
    }

    // Convert set to vector and then set bins for filter
    std::vector<int32_t> mats{mat_set.begin(), mat_set.end()};
    openmc_material_filter_set_bins(index_filter, mats.size(), mats.data());

    // Create tally and assign scores/filters
    openmc_extend_tallies(1, &openmcDriver.indexTally, nullptr);
    openmc_tally_set_type(openmcDriver.indexTally, "generic");
    openmc_tally_set_id(openmcDriver.indexTally, max_tally_id + 1);
    char score_array[][20]{"kappa-fission"};
    const char *scores[]{score_array[0]}; // OpenMC expects a const char**, ugh
    openmc_tally_set_scores(openmcDriver.indexTally, 1, scores);
    openmc_tally_set_filters(openmcDriver.indexTally, 1, &index_filter);
  }
}
