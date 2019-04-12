//! \file openmc_heat_driver.h
//! Driver for coupled surrogate heat transfer/OpenMC simulations
#ifndef ENRICO_OPENMC_HEAT_DRIVER_H
#define ENRICO_OPENMC_HEAT_DRIVER_H

#include "comm.h"
#include "coupled_driver.h"
#include "openmc_driver.h"
#include "heat_driver.h"

#include <mpi.h>

#include <unordered_map>
#include <vector>

namespace enrico {

class OpenmcHeatDriver : public CoupledDriver {

public:

  //! Initializes coupled OpenMC-surrogate driver with the given MPI communicator
  //!
  //! \param comm  The MPI communicator used for the coupled driver
  //! \param node  XML node containing settings
  explicit OpenmcHeatDriver(MPI_Comm comm, pugi::xml_node node);

  void update_heat_source() override;

  void update_temperature() override;

  Driver& getNeutronicsDriver() const override;

  Driver& getHeatDriver() const override;

  std::unique_ptr<OpenmcDriver> openmc_driver_; //! The OpenMC driver

  std::unique_ptr<SurrogateHeatDriver> heat_driver_; //! The heat surrogate driver

  // Mapping of surrogate rings to OpenMC cell instances and vice versa
  std::unordered_map<int, std::vector<int>> ring_to_cell_inst_;
  std::unordered_map<int, std::vector<int>> cell_inst_to_ring_;

private:

  //! Initialize mapping between OpenMC regions and surrogate fuel pin rings
  void init_mappings();

  //! Initialize tallies in OpenMC
  void init_tallies();
};

} // namespace enrico

#endif // ENRICO_OPENMC_HEAT_DRIVER_H
