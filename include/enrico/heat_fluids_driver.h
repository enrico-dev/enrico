//! \file heat_fluids_driver.h
//! Base class for single-physics heat-fluids solver
#ifndef HEAT_FLUIDS_DRIVER_H
#define HEAT_FLUIDS_DRIVER_H

#include "enrico/driver.h"
#include "enrico/geom.h"
#include "enrico/mpi_types.h"
#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

#include <cstddef> // for size_t

namespace enrico {

//! Base class for driver that controls a heat-fluids solve
class HeatFluidsDriver : public Driver {
public:
  HeatFluidsDriver(MPI_Comm comm, pugi::xml_node node);

  virtual ~HeatFluidsDriver() = default;

  //! Whether the calling rank has access to the fields returned
  //! by the 'temperature', 'density', and 'fluid_mask' methods
  //! \return Whether calling rank has access
  virtual bool has_coupling_data() const = 0;

  virtual int set_heat_source_at(int32_t local_elem, double heat) = 0;

  //! Return true if a local element is in the fluid region
  //! \param local_elem  A local element ID
  //! \return 1 if the local element is in fluid; 0 otherwise
  virtual int in_fluid_at(int32_t local_elem) const = 0;

  //! Get the number of local mesh elements
  //! \return Number of local mesh elements
  virtual int n_local_elem() const = 0;

  //! Get the number of global mesh elements
  //! \return Number of global mesh elements
  virtual std::size_t n_global_elem() const = 0;

  double pressure_bc_; //! System pressure in [MPa]

  //! Get temperature of local mesh elements
  //! \return Temperature of local mesh elements in [K]
  virtual std::vector<double> temperature() const = 0;

  //! Get density of local mesh elements
  //! \return Density of local mesh elements in [g/cm^3]
  virtual std::vector<double> density() const = 0;

  //! States whether each local region is in fluid
  //! \return For each local region, 1 if region is in fluid and 0 otherwise
  virtual std::vector<int> fluid_mask() const = 0;

  //! Get centroids of local mesh elements
  //! \return Centroids of local mesh elements
  virtual std::vector<Position> centroid() const = 0;

  //! Get volumes of local mesh elements
  //! \return Volumes of local mesh elements
  virtual std::vector<double> volume() const = 0;
};

} // namespace enrico

#endif // HEAT_FLUIDS_DRIVER_H
