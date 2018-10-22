#include "smrt/Assembly_Model.h"
#include "smrt/Shift_Solver.h"
#include "stream/message_passing.h"
#include "stream/nek_driver.h"
#include <mpi.h>

using AM = stream::Assembly_Model;

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  {
      std::string shift_filename = "singlerod_short.inp.xml";

      // Build assembly model
      std::vector<double> x_edges = {0.0, 1.26};
      std::vector<double> y_edges = {0.0, 1.26};
      std::vector<double> z_edges = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0};

      constexpr auto F = AM::FUEL;
      constexpr auto G = AM::GUIDE;

      std::vector<AM::PIN_TYPE> pins = {F};
      double height = 10.0;

      auto assembly = std::make_shared<AM>(pins, x_edges, y_edges, height);
      assembly->set_fuel_radius(0.418);
      assembly->set_clad_radius(0.475);
      assembly->set_guide_radius(0.61214);

      auto shift_solver = std::make_shared<stream::Shift_Solver>(
          assembly,
          shift_filename,
          z_edges);

      std::vector<double> power(5, 0.0);
      std::vector<double> fuel_temperature =
          {800.0, 1000.0, 1200.0, 1000.0, 800.0};
      std::vector<double> coolant_density =
          {0.74, 0.72, 0.70, 0.68, 0.66};

      shift_solver->solve(fuel_temperature, coolant_density, power);

      std::cout << "Shift computed power: ";
      for (const auto& val : power)
          std::cout << val << " ";
      std::cout << std::endl;

      stream::NekDriver nek_driver(MPI_COMM_WORLD);

      std::cout << "Done constructing Nek driver" << std::endl;

      std::cout << "Nek centroids: " << std::endl;
      for (int elem = 0; elem < nek_driver.nelgt_; ++elem) {
          auto c = nek_driver.get_global_elem_centroid(elem);
          std::cout << elem << " " << " " << c.x << " " << c.y
              << " " << c.z << std::endl;
      }

      if (nek_driver.active()) {
        nek_driver.init_step();
        nek_driver.solve_step();
        nek_driver.finalize_step();
      }
      nek_driver.comm_.Barrier();
  }

  MPI_Finalize();

  return 0;
}
