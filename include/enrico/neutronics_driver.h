//! \file neutronics_driver.h
//! Base class for single-physics neutronics solver
#ifndef NEUTRONICS_DRIVER_H
#define NEUTRONICS_DRIVER_H

#include "enrico/driver.h"
#include "enrico/geom.h"
#include "enrico/message_passing.h"

#include <gsl/gsl>
#include <xtensor/xtensor.hpp>

#include <vector>

namespace enrico {

using CellHandle = gsl::index;

//! Base class for driver that controls a neutronics solve
class NeutronicsDriver : public Driver {
public:
  explicit NeutronicsDriver(MPI_Comm comm)
    : Driver(comm)
  {}

  virtual ~NeutronicsDriver() = default;

  //! Get energy deposition in each material normalized to a given power
  //! \param power User-specified power in [W]
  //! \return Heat source in each material as [W/cm3]
  virtual xt::xtensor<double, 1> heat_source(double power) const = 0;

  //! Find cells corresponding to a vector of positions
  //! \param positions (x,y,z) coordinates to search for
  //! \return Handles to cells
  virtual std::vector<CellHandle> find(const std::vector<Position>& positions) = 0;

  //! Set the density of the material in a cell
  //! \param cell Handle to a cell
  //! \param rho Density in [g/cm^3]
  virtual void set_density(CellHandle cell, double rho) const = 0;

  //! Set the temperature of a cell
  //! \param cell Handle to a cell
  //! \param T Temperature in [K]
  virtual void set_temperature(CellHandle cell, double T) const = 0;

  //! Get the density of a cell
  //! \param cell Handle to a cell
  //! \return Cell density in [g/cm^3]
  virtual double get_density(CellHandle cell) const = 0;

  //! Get the temperature of a cell
  //! \param cell Handle to a cell
  //! \return Temperature in [K]
  virtual double get_temperature(CellHandle cell) const = 0;

  //! Get the volume of a cell
  //! \param cell Handle to a cell
  //! \return Volume in [cm^3]
  virtual double get_volume(CellHandle cell) const = 0;

  //! Detemrine whether a cell contains fissionable nuclides
  //! \param cell Handle to a cell
  //! \return Whether the cell contains fissionable nuclides
  virtual bool is_fissionable(CellHandle cell) const = 0;

  //! Determine number of cells participating in coupling
  //! \return Number of cells
  virtual std::size_t n_cells() const = 0;

  //! Create energy production tallies
  // TODO: Remove argument
  virtual void create_tallies(std::size_t n) = 0;
};

} // namespace enrico

#endif // NEUTRONICS_DRIVER_H
