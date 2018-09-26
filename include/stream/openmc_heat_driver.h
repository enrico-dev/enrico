//! \file openmc_heat_driver.h
//! Driver to run a surrogate finite difference heat equation solver
#ifndef STREAM_OPENMC_HEAT_DRIVER_H
#define STREAM_OPENMC_HEAT_DRIVER_H

#include <memory>
#include <unordered_map>
#include <vector>

#include <mpi.h>
#include <gsl/gsl>
#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

#include "base_drivers.h"
#include "openmc_driver.h"
#include "comm.h"
#include "heat_driver.h"
#include "geom.h"

namespace stream {

class SurrogateHeatDriver : public HeatFluidsDriver {
public:
  explicit SurrogateHeatDriver(MPI_Comm comm, pugi::xml_node node);
  ~SurrogateHeatDriver() { };

  void init_step() { };
  void solve_step() { };
  void finalize_step() { };

  void set_heat_source();

  std::unique_ptr<HeatSolver> solver_;

  // Data on fuel pins
  xt::xtensor<double, 2> pin_centers_; //!< (x,y) values for center of fuel pins
  xt::xtensor<double, 1> z_; //!< Bounding z-values for axial segments
};

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
  double power_;

  // Mapping of surrogate rings to OpenMC cell instances and vice versa
  std::unordered_map<int, std::vector<int>> ring_to_cell_inst_;
  std::unordered_map<int, std::vector<int>> cell_inst_to_ring_;

private:
  void init_mappings();
  void init_tallies();
};

} // namespace stream

#endif // STREAM_OPENMC_HEAT_DRIVER_H
