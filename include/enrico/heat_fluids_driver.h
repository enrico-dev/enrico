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

  //! Get the temperature in each region
  //! \return Temperature in each region as [K]
  xt::xtensor<double, 1> temperature() const;

  //! Get the density in each region. The interpretation of this density,
  //! i.e. whether it refers to fluid elements or both fluid and solid
  //! elements, is to the discretion of the particular driver.
  //! \return Density in each region as [g/cm^3]
  xt::xtensor<double, 1> density() const;

  //! States whether each region is in fluid
  //! \return For each region, 1 if region is in fluid and 0 otherwise
  std::vector<int> fluid_mask() const;

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

  //! Get the centroids of all mesh elements
  //! \return Vector of all centroids
  std::vector<Position> centroids() const;

  //! Get the volumes of all mesh elements
  //! \return Vector of all volumes
  std::vector<double> volumes() const;

  double pressure_bc_; //! System pressure in [MPa]

  //! The displacements of local elements, relative to rank 0. Used in an MPI
  //! Gatherv operation.
  // TODO: Move to private
  std::vector<int32_t> local_displs_;

  //! The number of local elements in each rank.
  // TODO: Move to private
  std::vector<int32_t> local_counts_;

  //! Get temperature of local mesh elements
  //! \return Temperature of local mesh elements in [K]
  virtual std::vector<double> temperature_local() const = 0;

  //! Get density of local mesh elements
  //! \return Density of local mesh elements in [g/cm^3]
  virtual std::vector<double> density_local() const = 0;

  //! States whether each local region is in fluid
  //! \return For each local region, 1 if region is in fluid and 0 otherwise
  virtual std::vector<int> fluid_mask_local() const = 0;

  //! Get centroids of local mesh elements
  //! \return Centroids of local mesh elements
  virtual std::vector<Position> centroid_local() const = 0;

  //! Get volumes of local mesh elements
  //! \return Volumes of local mesh elements
  virtual std::vector<double> volume_local() const = 0;

protected:
  //! Initialize the counts and displacements of local elements for each MPI Rank.
  void init_displs();

private:
  //! Gather local distributed field into global field (on rank 0)
  //! \return Global field collected from all ranks
  template<typename T>
  std::vector<T> gather(const std::vector<T>& local_field) const;
};

template<typename T>
std::vector<T> HeatFluidsDriver::gather(const std::vector<T>& local_field) const
{
  std::vector<T> global_field;

  if (this->active()) {
    if (this->has_coupling_data()) {
      global_field.resize(this->n_global_elem());
    }

    // Gather all the local quantities on to the root process
    comm_.Gatherv(local_field.data(),
                  local_field.size(),
                  get_mpi_type<T>(),
                  global_field.data(),
                  local_counts_.data(),
                  local_displs_.data(),
                  get_mpi_type<T>());
  }

  return global_field;
}

} // namespace enrico

#endif // HEAT_FLUIDS_DRIVER_H
