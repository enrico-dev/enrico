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

#include <unordered_map>
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

  //////////////////////////////////////////////////////////////////////////////
  // NeutronicsDriver interface

  //! Find cells corresponding to a vector of positions
  //! \param positions (x,y,z) coordinates to search for
  //! \return Handles to cells
  std::vector<CellHandle> find(const std::vector<Position>& position) override;

  //! Set the density of the material in a cell
  //! \param cell Handle to a cell
  //! \param rho Density in [g/cm^3]
  void set_density(CellHandle cell, double rho) const override;

  //! Set the temperature of a cell
  //! \param cell Handle to a cell
  //! \param T Temperature in [K]
  void set_temperature(CellHandle cell, double T) const override;

  //! Get the density of a cell
  //! \param cell Handle to a cell
  //! \return Cell density in [g/cm^3]
  double get_density(CellHandle cell) const override;

  //! Get the temperature of a cell
  //! \param cell Handle to a cell
  //! \return Temperature in [K]
  double get_temperature(CellHandle cell) const override;

  //! Get the volume of a cell
  //! \param cell Handle to a cell
  //! \return Volume in [cm^3]
  double get_volume(CellHandle cell) const override;

  //! Detemrine whether a cell contains fissionable nuclides
  //! \param cell Handle to a cell
  //! \return Whether the cell contains fissionable nuclides
  bool is_fissionable(CellHandle cell) const override;

  std::size_t n_cells() const override { return cells_.size(); }

  //! Create energy production tallies
  void create_tallies() override;

  //! Determine number of cells participating in coupling
  //! \return Number of cells
  xt::xtensor<double, 1> heat_source(double power) const final;

  std::string cell_label(CellHandle cell) const;

  gsl::index cell_index(CellHandle cell) const override;

  //////////////////////////////////////////////////////////////////////////////
  // Driver interface

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

private:
  CellInstance& cell_instance(CellHandle cell);
  const CellInstance& cell_instance(CellHandle cell) const;

  // Data members
  openmc::Tally* tally_;               //!< Fission energy deposition tally
  openmc::CellInstanceFilter* filter_; //!< Cell instance filter
  std::vector<CellInstance> cells_;    //!< Array of cell instances
  std::unordered_map<CellHandle, gsl::index>
    cell_index_;            //!< Map handles to index in cells_
  int n_fissionable_cells_; //!< Number of fissionable cells in model
};

} // namespace enrico

#endif // ENRICO_OPENMC_DRIVER_H
