//! \file openmc_interface.h
//! Classes to access OpenMC data
#ifndef STREAM_OPENMC_INTERFACE_H
#define STREAM_OPENMC_INTERFACE_H

#include <cstdint>

#include "openmc/material.h"
#include "stream_geom.h"

namespace stream {

//! Get/set a cell's data, including data linked to its material
class CellInstance {
public:

  //! Given a position, find the cell and material IDs of the bounding cell
  //!
  //! The position should be in OpenMC's units (cm).
  //!
  //! \param position The coordinate for the desired cell
  explicit CellInstance(Position position);

  //! Get this cell's material
  //! \return A pointer to the Material associated with this cell
  openmc::Material* material() const;

  //! Set the density of this cell's materials
  void set_density(double rho) const;

  //! Set the temperature of this cell
  void set_temperature(double T) const;

  int32_t index_; //!< Index in global cells array
  int32_t instance_; //!< Index of cell instance
  int32_t material_index_; //!< Index of material in this instance
  double volume_ {0.0}; //!< volume of cell instance in [cm^3]
};

}

#endif // STREAM_OPENMC_INTERFACE_H
