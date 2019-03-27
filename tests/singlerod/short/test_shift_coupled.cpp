#include <string>
#include <vector>

#include "Nemesis/comm/global.hh"
#include "smrt/Assembly_Model.h"
#include "smrt/Coupled_Solver.h"
#include <mpi.h>

using AM = enrico::Assembly_Model;

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  {
  std::string shift_filename  = "singlerod_short.inp.xml";
  std::string enrico_filename = "enrico.xml";

  // Build assembly model
  std::vector<double> x_edges = {-0.63, 0.63};
  std::vector<double> y_edges = {-0.63, 0.63};
  std::vector<double> z_edges = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0};

  constexpr auto F = AM::FUEL;
  constexpr auto G = AM::GUIDE;

  std::vector<AM::PIN_TYPE> pins = {F};
  double height = z_edges.back();

  auto assembly = std::make_shared<AM>(pins, x_edges, y_edges, height);
  assembly->set_fuel_radius(0.418);
  assembly->set_clad_radius(0.475);
  assembly->set_guide_radius(0.61214);

  double power_norm = 1800.0;
  auto coupled_solver = std::make_shared<enrico::Coupled_Solver>(
      assembly,
      z_edges,
      shift_filename,
      enrico_filename,
      power_norm,
      MPI_COMM_WORLD,
      MPI_COMM_WORLD);

  coupled_solver->solve();
  }

  MPI_Finalize();

  return 0;
}
