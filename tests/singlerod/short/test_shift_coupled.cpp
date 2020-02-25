#include <string>
#include <vector>

#include "Nemesis/comm/global.hh"
#include "enrico/mpi_types.h"
#include "smrt/Assembly_Model.h"
#include "smrt/shift_nek_driver.h"
#include <mpi.h>

using AM = enrico::Assembly_Model;

int main(int argc, char* argv[])
{

  MPI_Init(&argc, &argv);
  enrico::init_mpi_datatypes();

  {
    std::string shift_filename = "singlerod_short.inp.xml";

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

    auto coupled_solver = std::make_shared<enrico::ShiftNekDriver>(
      assembly, z_edges, shift_filename, MPI_COMM_WORLD, MPI_COMM_WORLD);

    coupled_solver->solve();
  }

  enrico::free_mpi_datatypes();
  MPI_Finalize();

  return 0;
}
