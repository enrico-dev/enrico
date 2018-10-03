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
  explicit OpenmcHeatDriver(MPI_Comm comm, pugi::xml_node node);
  ~OpenmcHeatDriver() { };

  void solve_step();

  void update_heat_source();

  void update_temperature();

  // Data
  Comm comm_;
  Comm intranode_comm_;
  std::unique_ptr<OpenmcDriver> openmc_driver_;
  std::unique_ptr<SurrogateHeatDriver> heat_driver_;
  double power_; //!< Power in [W]

  // Mapping of surrogate rings to OpenMC cell instances and vice versa
  std::unordered_map<int, std::vector<int>> ring_to_cell_inst_;
  std::unordered_map<int, std::vector<int>> cell_inst_to_ring_;

private:
  void init_mappings();
  void init_tallies();
};

} // namespace stream

#endif // STREAM_OPENMC_HEAT_DRIVER_H
