//! \file openmc_heat_driver.h
//! Driver to run a surrogate finite difference heat equation solver
#ifndef STREAM_OPENMC_HEAT_DRIVER_H
#define STREAM_OPENMC_HEAT_DRIVER_H

#include "comm.h"
#include "openmc_driver.h"
#include "heat_driver.h"

#include <mpi.h>
#include "pugixml.hpp"

#include <memory>
#include <unordered_map>
#include <vector>

namespace stream {

class OpenmcHeatDriver {
public:
  //! Initializes coupled OpenMC-surrogate driver with the given MPI communicator
  //!
  //! \param comm  The MPI communicator used for the coupled driver
  //! \param node  XML node containing settings
  explicit OpenmcHeatDriver(MPI_Comm comm, pugi::xml_node node);

  //! Destroy the coupled driver
  ~OpenmcHeatDriver() { };

  //! Solve a single coupled step
  //!
  //! This method runs OpenMC to obtain a heat source, updates the heat source
  //! for the surrogate heat-fluids driver, solves the heat equation, and then
  //! updates the temperature for cell instances in OpenMC.
  void solve_step();

  //! Update the heat source for the surrogate heat-fluids driver
  void update_heat_source();

  //! Update the temperature in each region for OpenMC
  void update_temperature();

  //! Run one timstep
  void solve_in_time();

  //! Write driver data to VTK file
  void to_vtk(int iteration = -1, int timestep = -1) {
    heat_driver_->to_vtk(iteration, timestep);
  }

  // Data
  Comm comm_;  //!< The communicator used to run this driver
  std::unique_ptr<OpenmcDriver> openmc_driver_; //! The OpenMC driver
  std::unique_ptr<SurrogateHeatDriver> heat_driver_; //! The heat surrogate driver
  double power_; //!< Power in [W]
  int max_timesteps_; //! Maximum of timesteps
  int max_picard_iter_; //! Maximum number of Picard iterations per timestep

  // Mapping of surrogate rings to OpenMC cell instances and vice versa
  std::unordered_map<int, std::vector<int>> ring_to_cell_inst_;
  std::unordered_map<int, std::vector<int>> cell_inst_to_ring_;

private:
  //! Initialize mapping between OpenMC regions and surrogate fuel pin rings
  void init_mappings();

  //! Initialize tallies in OpenMC
  void init_tallies();
};

} // namespace stream

#endif // STREAM_OPENMC_HEAT_DRIVER_H
