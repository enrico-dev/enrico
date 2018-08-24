#ifndef STREAM_OPENMC_NEK_DRIVER_H
#define STREAM_OPENMC_NEK_DRIVER_H

#include "drivers.h"
#include "mpi.h"
#include "nek_interface.h"
#include "nek_driver.h"
#include "openmc/capi.h"
#include "stream_geom.h"
#include "stream_const.h"

#include <algorithm> // for max, fill, copy
#include <iterator> // for back_inserter
#include <unordered_set>

namespace stream {

// This is not actually derived from CoupledDriver.  Currently, it is unclear
// how or if the base class will be implemented.  The issue will be revisited
class OpenmcNekDriver {
public:
  OpenmcNekDriver(int argc, char *argv[], MPI_Comm coupled_comm, MPI_Comm openmc_comm,
                  MPI_Comm nek_comm, MPI_Comm intranode_comm);
  ~OpenmcNekDriver() {};

  void update_heat_source();
  void update_temperature();

  Comm comm_;
  Comm intranode_comm_;
  OpenmcDriver openmc_driver_;
  NekDriver nek_driver_;
private:
  void init_mappings();
  void init_tallies();

  int get_heat_index(int32_t mat_index) const {
    return heat_index_.at(mat_index - 1);
  }

  // Map that gives a list of Nek element global indices for a given OpenMC
  // material index
  std::unordered_map<int32_t, std::vector<int>> mat_to_elems_;

  // Map that gives the OpenMC material index for a given Nek global element index
  std::unordered_map<int, int32_t> elem_to_mat_;

  // Mapping of material indices (minus 1) to positions in array of heat sources that
  // is used during update_heat_source
  std::vector<int> heat_index_;

  int32_t n_materials_;
};

} // namespace stream

#endif //STREAM_OPENMC_NEK_DRIVER_H
