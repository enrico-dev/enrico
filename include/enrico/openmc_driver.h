//! \file openmc_driver.h
//! Driver to initialize and run OpenMC in stages
#ifndef ENRICO_OPENMC_DRIVER_H
#define ENRICO_OPENMC_DRIVER_H

#include "driver.h"
#include "openmc_interface.h"
#include "geom.h"

#include <gsl/gsl>
#include <mpi.h>
#include "xtensor/xtensor.hpp"
#include "openmc/tallies/tally.h"

#include <vector>

namespace enrico {

//! Driver to initialize and run OpenMC in stages
class OpenmcDriver : public Driver {
public:
  //! One-time initalization of OpenMC and member variables
  //! \param comm An existing MPI communicator used to inialize OpenMC
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

  //! Initialization required in each Picard iteration
  virtual void init_step() override;

  //! Runs OpenMC for one Picard iteration
  virtual void solve_step() override;

  //! Writes OpenMC output for given timestep and iteration
  //! \param timestep timestep index
  //! \param iteration iteration index
  virtual void write_step(int timestep = -1, int iteration = -1) override;

  //! Finalization required in each Picard iteration
  virtual void finalize_step() override;

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

} // namespace enrico

#endif //ENRICO_OPENMC_DRIVER_H
