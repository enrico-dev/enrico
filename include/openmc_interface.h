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

  int32_t index_;
  int32_t instance_;
  int32_t material_index_;
};


}

#endif // STREAM_OPENMC_INTERFACE_H
