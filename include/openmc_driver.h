//! \file openmc_driver.h
//! Driver to initialize and run OpenMC in stages
#ifndef STREAM_OPENMC_DRIVER_H
#define STREAM_OPENMC_DRIVER_H

#include "base_drivers.h"
#include "mpi.h"
#include "openmc_interface.h"
#include "stream_geom.h"
#include <vector>

namespace stream {

//! Driver to initialize and run OpenMC in stages
class OpenmcDriver : public TransportDriver {
public:
  //! One-time initalization of OpenMC and member variables
  //! \param argc Number of command-line arguments
  //! \param argv Values of command-line arguments
  //! \param comm An exiting MPI communicator used to inialize OpenMC
  OpenmcDriver(int argc, char* argv[], MPI_Comm comm);

  //! One-time finalization of OpenMC
  ~OpenmcDriver();

  //! Initalization required in each Picard iteration
  void init_step();

  //! Runs OpenMC for one Picard iteration
  void solve_step();

  //! Finaliztion required in each Picard iteration
  void finalize_step();

  //! Get the coordinate of a given material's centroid
  //!
  //! This coordinate has the units used by OpenMC (cm)
  //!
  //! \param mat_id A material ID
  //! \return The coordinate of the material's centroid in cm
  Position get_mat_centroid(int32_t mat_id) const;

  // Data
  int32_t index_tally_;   //!< Index in tallies array for fission tally
  int32_t index_filter_;  //!< Index in filters arrays for material filter
  std::vector<CellInstance> cells_;  //!< Array of cell instances
};

} // namespace stream

#endif //STREAM_OPENMC_DRIVER_H
