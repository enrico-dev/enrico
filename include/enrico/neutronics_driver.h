//! \file neutronics_driver.h
//! Base class for single-physics neutronics solver
#ifndef NEUTRONICS_DRIVER_H
#define NEUTRONICS_DRIVER_H

#include "enrico/cell_handle.h"
#include "enrico/driver.h"
#include "enrico/geom.h"
#include "enrico/mpi_types.h"

#include <gsl/gsl>
#include <xtensor/xtensor.hpp>

#include <vector>

namespace enrico {

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

  //! Get the k-effective of a run
  virtual double get_k_effective() const = 0;

  //! Get the boron concentration
  virtual double get_boron_ppm() const = 0;

  //! Get the Boronated H2O density
  virtual double get_H2O_dens() const = 0;

  //! Set the Boron concentration in a cell
  //! \param ppm Boric acid concentration in [ppm] !TODO: by wgt?
  //! \param H2Odens water density in [g/cm^3]
  virtual double set_boron_ppm(double ppm, double H2Odens) const = 0;

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
  virtual void create_tallies() = 0;

  //! Get a label for a cell
  //! \param cell Handle to a clel
  //! \return Label for the cell
  virtual std::string cell_label(CellHandle cell) const = 0;

  //! Get the index of the given handle in the cells_ ordered mapping.
  //!
  //! Currently, this is used to infer the index in the heat source array that corresponds
  //! to a given cell handle.
  //!
  //! \param cell An existing cell handle
  //! \return The index of the handle in the cells_ ordered mapping
  virtual gsl::index cell_index(CellHandle cell) const = 0;
};

} // namespace enrico

#endif // NEUTRONICS_DRIVER_H
