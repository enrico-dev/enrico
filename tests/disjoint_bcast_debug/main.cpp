#include "enrico/comm.h"
#include "enrico/comm_split.h"
#include "mpi.h"
#include "pugixml.hpp"

#include <array>

class TestDriver {
public:
  enrico::Comm left_comm, right_comm;

  TestDriver()
  {
    std::array<int, 2> nodes{2, 2};
    std::array<int, 2> procs_per_node{4, 1};

    enrico::Comm super_comm(MPI_COMM_WORLD);
    std::array<enrico::Comm, 2> driver_comms;
    enrico::Comm intranode_comm, coupling_comm;

    enrico::get_driver_comms(
      super_comm, nodes, procs_per_node, driver_comms, intranode_comm, coupling_comm);

    left_comm = driver_comms[0];
    right_comm = driver_comms[1];
  }

  void test_right_bcast()
  {
    if (right_comm.comm != MPI_COMM_NULL) {
      double val = right_comm.rank == 0 ? 3.14 : 0.0;
      MPI_Bcast(&val, 1, MPI_DOUBLE, 0, right_comm.comm);

      right_comm.message("Rank " + std::to_string(right_comm.rank) +
                         ": val = " + std::to_string(val));
    }
  }
};

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  TestDriver d;
  d.test_right_bcast();

  MPI_Finalize();
}
