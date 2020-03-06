#include "enrico/comm.h"
#include "enrico/comm_split.h"
#include "pugixml.hpp"

#include <mpi.h>
#include <unistd.h>

class DisjointDriver {
public:
  DisjointDriver(MPI_Comm comm, pugi::xml_node node);

  enrico::Comm super_comm;

  std::array<enrico::Comm, 2> driver_comms;
  enrico::Comm intranode_comm, coupling_comm;
};

DisjointDriver::DisjointDriver(MPI_Comm comm, pugi::xml_node node)
  : super_comm(comm)
{
  std::array<int, 2> nodes{node.child("neutronics_nodes").text().as_int(),
                           node.child("heat_fluids_nodes").text().as_int()};
  std::array<int, 2> procs_per_node{
    node.child("neutronics_procs_per_node").text().as_int(),
    node.child("heat_fluids_procs_per_node").text().as_int()};
  enrico::get_driver_comms(
    super_comm, nodes, procs_per_node, driver_comms, intranode_comm, coupling_comm);
}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  char hostname[_POSIX_HOST_NAME_MAX];
  gethostname(hostname, _POSIX_HOST_NAME_MAX);

  // Parse enrico.xml file
  pugi::xml_document doc;
  auto result = doc.load_file("enrico.xml");
  if (!result) {
    throw std::runtime_error{"Unable to load enrico.xml file"};
  }
  // Get root element
  auto root = doc.document_element();

  enrico::Comm world(MPI_COMM_WORLD);

  DisjointDriver driver(world.comm, root);
  for (int i = 0; i < world.size; ++i) {
    if (world.rank == i) {
      if (i == 0) {
        std::cout << "Host\t\tWorld\tLeft\tRight\tIntra\tCoup" << std::endl;
      }
      std::cout << hostname << "\t" << world.rank << "\t" << driver.driver_comms[0].rank
                << "\t" << driver.driver_comms[1].rank << "\t"
                << driver.intranode_comm.rank << "\t" << driver.coupling_comm.rank << "\t"
                << std::endl;
    }
    MPI_Barrier(world.comm);
  }

  MPI_Finalize();
}