//! \file openmc_interface.h
//! Class to get/set a cell's data, including data linked to its material
#ifndef ENRICO_CELL_INSTANCE_H
#define ENRICO_CELL_INSTANCE_H

#include <cstdint>
#include <functional> // for hash

#include "geom.h"
#include "openmc/material.h"

namespace enrico {

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

  //! Get the temperature of this cell
  //! \return cell temperature in [K]
  double get_temperature() const;

  //! Get the density of this cell
  //! \return cell density in [g/cm^3]
  double get_density() const;

  //! Check for equality
  bool operator==(const CellInstance& other) const;

  int32_t index_;          //!< Index in global cells array
  int32_t instance_;       //!< Index of cell instance
  int32_t material_index_; //!< Index of material in this instance
  double volume_{0.0};     //!< volume of cell instance in [cm^3]
};

} // namespace enrico

namespace std {

template<>
struct hash<enrico::CellInstance> {
  // Taken from https://stackoverflow.com/a/17017281
  std::size_t operator()(const enrico::CellInstance& k) const
  {
    std::size_t res = 17;
    res = 31 * res + std::hash<int32_t>()(k.index_);
    res = 31 * res + std::hash<int32_t>()(k.instance_);
    return res;
  }
};

} // namespace std

#endif // ENRICO_CELL_INSTANCE_H
