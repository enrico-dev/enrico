#include "enrico/comm.h"
#include "enrico/comm_split.h"
#include "mpi.h"
#include "pugixml.hpp"

#include <array>
#include <iostream>
#include <iomanip>
#include <unistd.h>

class TestDriver {
public:
  enrico::Comm left, right;

  TestDriver()
  {
    std::array<int, 2> nodes{2, 2};
    std::array<int, 2> procs_per_node{1, 4};

    enrico::Comm super_comm(MPI_COMM_WORLD);
    std::array<enrico::Comm, 2> driver_comms;
    enrico::Comm intranode_comm, coupling_comm;

    enrico::get_driver_comms(
      super_comm, nodes, procs_per_node, driver_comms, intranode_comm, coupling_comm);

    left = driver_comms[0];
    right = driver_comms[1];
  }

  void test_bcast(enrico::Comm c)
  {
    if (c.comm != MPI_COMM_NULL) {
      double val = c.rank == 0 ? 3.14 : 0.0;
      MPI_Bcast(&val, 1, MPI_DOUBLE, 0, c.comm);

      for (int i = 0; i < c.size; ++i) {
        c.message("Rank " + std::to_string(c.rank) + ": val = " + std::to_string(val),
                           i);
      }
    }
  }

  void comm_report()
  {
    char c[_POSIX_HOST_NAME_MAX];
    gethostname(c, _POSIX_HOST_NAME_MAX);
    std::string hostname{c};

    enrico::Comm world(MPI_COMM_WORLD);

    // Padding for fields
    int hostw = std::max(8UL, hostname.size()) + 2;
    int rankw = 7;

    for (int i = 0; i < world.size; ++i) {
      if (world.rank == i) {
        if (i == 0) {
          std::cout << std::left << std::setw(hostw) << "Hostname" << std::right
                    << std::setw(rankw) << "World" << std::right 
                    << std::setw(rankw) << "Left" << std::right 
                    << std::setw(rankw) << "Right" << std::endl;
        }
        std::cout << std::left << std::setw(hostw) << hostname << std::right
                  << std::setw(rankw) << world.rank << std::right 
                  << std::setw(rankw) << left.rank << std::right 
                  << std::setw(rankw) << right.rank << std::endl;
      }
      MPI_Barrier(world.comm);
    }
  }
};

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  TestDriver d;
  d.comm_report();
  d.test_bcast(d.left);

  MPI_Finalize();
}
