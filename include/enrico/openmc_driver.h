//! \file openmc_driver.h
//! Driver to initialize and run OpenMC in stages
#ifndef ENRICO_OPENMC_DRIVER_H
#define ENRICO_OPENMC_DRIVER_H

#include "cell_instance.h"
#include "geom.h"
#include "neutronics_driver.h"

#include "openmc/tallies/tally.h"
#include <gsl/gsl>
#include <mpi.h>

#include <vector>

namespace enrico {

//! Driver to initialize and run OpenMC in stages
class OpenmcDriver : public NeutronicsDriver {
public:
  //! One-time initalization of OpenMC and member variables
  //! \param comm An existing MPI communicator used to inialize OpenMC
  explicit OpenmcDriver(MPI_Comm comm);

  //! One-time finalization of OpenMC
  ~OpenmcDriver();

  //! Create energy production tallies for a list of materials
  //! \param[in] materials  Indices into OpenMC's materials array
  void create_tallies(gsl::span<int32_t> materials);

  xt::xtensor<double, 1> heat_source(double power) const final;

  //! Initialization required in each Picard iteration
  void init_step() final;

  //! Runs OpenMC for one Picard iteration
  void solve_step() final;

  //! Writes OpenMC output for given timestep and iteration
  //! \param timestep timestep index
  //! \param iteration iteration index
  void write_step(int timestep, int iteration) final;

  //! Finalization required in each Picard iteration
  void finalize_step() final;

  // Data
  openmc::Tally* tally_;            //!< Fission energy deposition tally
  int32_t index_filter_;            //!< Index in filters arrays for material filter
  std::vector<CellInstance> cells_; //!< Array of cell instances
};

} // namespace enrico

#endif // ENRICO_OPENMC_DRIVER_H
