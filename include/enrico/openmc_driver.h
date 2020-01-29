//! \file openmc_driver.h
//! Driver to initialize and run OpenMC in stages
#ifndef ENRICO_OPENMC_DRIVER_H
#define ENRICO_OPENMC_DRIVER_H

#include "enrico/cell_instance.h"
#include "enrico/geom.h"
#include "enrico/neutronics_driver.h"

#include "openmc/cell.h"
#include "openmc/tallies/filter_cell_instance.h"
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

  std::vector<CellHandle> find(const std::vector<Position>& position) override;
  void set_density(CellHandle cell, double rho) const override;
  void set_temperature(CellHandle cell, double T) const override;
  double get_density(CellHandle cell) const override;
  double get_temperature(CellHandle cell) const override;
  double get_volume(CellHandle cell) const override;
  bool is_fissionable(CellHandle cell) const override;
  std::size_t n_cells() const override { return cells_.size(); }

  //! Create energy production tallies for a list of cell instances
  //! \param[in] cells  Sequence of OpenMC cell instances
  void create_tallies(std::size_t n) override;

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
  openmc::Tally* tally_;               //!< Fission energy deposition tally
  openmc::CellInstanceFilter* filter_; //!< Cell instance filter
  std::vector<CellInstance> cells_;    //!< Array of cell instances
  int n_fissionable_cells_;            //!< Number of fissionable cells in model
};

} // namespace enrico

#endif // ENRICO_OPENMC_DRIVER_H
