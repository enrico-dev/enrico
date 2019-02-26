//! \file openmc_driver.h
//! Driver to initialize and run OpenMC in stages
#ifndef STREAM_OPENMC_DRIVER_H
#define STREAM_OPENMC_DRIVER_H

#include "base_drivers.h"
#include "openmc_interface.h"
#include "geom.h"

#include <gsl/gsl>
#include <mpi.h>
#include "xtensor/xtensor.hpp"
#include "openmc/tallies/tally.h"

#include <vector>

namespace stream {

//! Driver to initialize and run OpenMC in stages
class OpenmcDriver : public TransportDriver {
public:
  //! One-time initalization of OpenMC and member variables
  //! \param comm An exiting MPI communicator used to inialize OpenMC
  OpenmcDriver(MPI_Comm comm);

  //! One-time finalization of OpenMC
  ~OpenmcDriver();

  //! Create energy production tallies for a list of materials
  //! \param[in] materials  Indices into OpenMC's materials array
  void create_tallies(gsl::span<int32_t> materials);

  //! Get energy deposition in each material
  //!
  //! \param power User-specified power in [W]
  //! \return Heat source in each material as [W/cm3]
  xt::xtensor<double, 1> heat_source(double power);

  //! Initalization required in each Picard iteration
  void init_step();

  //! Runs OpenMC for one Picard iteration
  //!
  //! \param i Iteration number
  void solve_step(int i);

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
  openmc::Tally* tally_; //!< Fission energy deposition tally
  int32_t index_filter_;  //!< Index in filters arrays for material filter
  std::vector<CellInstance> cells_;  //!< Array of cell instances
};

} // namespace stream

#endif //STREAM_OPENMC_DRIVER_H
