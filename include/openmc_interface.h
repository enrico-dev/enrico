#ifndef STREAM_OPENMC_INTERFACE_H
#define STREAM_OPENMC_INTERFACE_H

#include <cstdint>

#include "openmc/material.h"
#include "stream_geom.h"

namespace stream {

class CellInstance {
public:
  explicit CellInstance(Position position);

  openmc::Material* material() const;
  void set_density(double rho) const;
  void set_temperature(double T) const;

  int32_t index_; //!< Index in global cells array
  int32_t instance_; //!< Index of cell instance
  int32_t material_index_; //!< Index of material in this instance
  double volume_; //!< volume of cell instance in [cm^3]
};


}

#endif // STREAM_OPENMC_INTERFACE_H
