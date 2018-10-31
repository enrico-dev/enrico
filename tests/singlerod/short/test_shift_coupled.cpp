#include "Nemesis/comm/global.hh"
#include "smrt/Assembly_Model.h"
#include "smrt/Shift_Solver.h"
#include "stream/message_passing.h"
#include "stream/nek_driver.h"
#include "stream/nek_interface.h"
#include "stream/error.h"
#include <mpi.h>

using AM = stream::Assembly_Model;

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  {
      std::string shift_filename = "singlerod_short.inp.xml";

      // Build assembly model
      std::vector<double> x_edges = {-0.63, 0.63};
      std::vector<double> y_edges = {-0.63, 0.63};
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

      stream::NekDriver nek_driver(MPI_COMM_WORLD);

      if (nek_driver.active()) {
        nek_driver.init_step();
        nek_driver.solve_step();
        nek_driver.finalize_step();
        nek_driver.comm_.Barrier();

        int num_local  = nek_driver.nelt_;
        int num_global = nek_driver.nelgt_;

        std::vector<stream::Position> local_centroids(num_local);
        std::vector<double> local_volumes(num_local);
        std::vector<double> local_temperatures(num_local);
        for (int elem = 0; elem < num_local; ++elem)
        {
            // Need to give Nek the Fortran (1-based) index
            local_centroids[elem] = nek_driver.get_local_elem_centroid(elem+1);
            local_volumes[elem]   = nek_driver.get_local_elem_volume(elem+1);
            local_temperatures[elem] = nek_driver.get_local_elem_temperature(elem+1);
        }

        // set up MPI datatype
        // Currently, this sets up only position_mpi_datatype
        MPI_Datatype position_mpi_datatype;
        {
            stream::Position p;
            int blockcounts[3] = {1, 1, 1};
            MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
            MPI_Aint displs[3];

            // Get displacements of struct members
            MPI_Get_address(&p.x, &displs[0]);
            MPI_Get_address(&p.y, &displs[1]);
            MPI_Get_address(&p.z, &displs[2]);

            // Make the displacements relative
            displs[2] -= displs[0];
            displs[1] -= displs[0];
            displs[0] = 0;

            // Make datatype
            MPI_Type_create_struct(3, blockcounts, displs, types, &position_mpi_datatype);
            MPI_Type_commit(&position_mpi_datatype);
        }

        // Gather global centroid list to root process
        std::vector<stream::Position> global_centroids(num_global);
        nek_driver.comm_.Gatherv(local_centroids.data(),
                                 num_local,
                                 position_mpi_datatype,
                                 global_centroids.data(),
                                 nek_driver.local_counts_.data(),
                                 nek_driver.local_displs_.data(),
                                 position_mpi_datatype);

        // Broadcast centroids
        nek_driver.comm_.Bcast(global_centroids.data(),
                               num_global,
                               position_mpi_datatype);

        // Repeat for volumes
        std::vector<double> global_volumes(num_global);
        nek_driver.comm_.Gatherv(local_volumes.data(),
                                 num_local,
                                 MPI_DOUBLE,
                                 global_volumes.data(),
                                 nek_driver.local_counts_.data(),
                                 nek_driver.local_displs_.data(),
                                 MPI_DOUBLE);

        // Broadcast centroids
        nek_driver.comm_.Bcast(global_volumes.data(),
                               num_global,
                               MPI_DOUBLE);

        // Register centroids with Shift solver
        shift_solver->set_centroids_and_volumes(
            global_centroids, global_volumes);

        std::vector<double> global_temperature(num_global, 1000.0);
        std::vector<double> global_density(num_global, 0.75);
        std::vector<double> global_power(num_global, 0.0);

        shift_solver->solve(global_temperature, global_density, global_power);

        for (auto& val : global_power)
            val *= 1600.0;

        // Update heat source in Nek and solve again
        int displacement = nek_driver.local_displs_[nemesis::node()];
        for (int elem = 0; elem < num_local; ++elem)
        {
            stream::err_chk(
                nek_set_heat_source(elem, global_power[elem+displacement]),
                    "Error setting heat source for local element " +
                        std::to_string(elem));
        }

        nek_driver.solve_step();
      }
  }

  MPI_Finalize();

  return 0;
}
