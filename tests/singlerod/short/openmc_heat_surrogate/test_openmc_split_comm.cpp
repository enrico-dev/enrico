#include "enrico/comm_split.h"
#include "enrico/openmc_driver.h"
#include "enrico/surrogate_heat_driver.h"
#include <iomanip>
#include <mpi.h>
#include <unistd.h>

namespace enrico {

class TestDriver {
public:

  TestDriver(MPI_Comm comm, pugi::xml_node node)
    : comm_(comm)
  {
    auto neut_node = node.child("neutronics");
    auto heat_node = node.child("heat_fluids");
    auto coup_node = node.child("coupling");

    // Create communicators
    std::array<int, 2> nodes{neut_node.child("nodes").text().as_int(),
                             heat_node.child("nodes").text().as_int()};
    std::array<int, 2> procs_per_node{neut_node.child("procs_per_node").text().as_int(),
                                      heat_node.child("procs_per_node").text().as_int()};
    std::array<Comm, 2> driver_comms;
    Comm intranode_comm; // Not used in current comm scheme
    Comm coupling_comm;  // Not used in current comm scheme

    get_driver_comms(
      comm_, nodes, procs_per_node, driver_comms, intranode_comm, coupling_comm);

    auto neutronics_comm = driver_comms[0];
    auto heat_comm = driver_comms[1];

    // Instantiate neutronics driver
    std::string neut_driver = neut_node.child_value("driver");
    neutronics_driver_ = std::make_unique<OpenmcDriver>(neutronics_comm.comm);

    // Instantiate heat-fluids driver
    std::string s = heat_node.child_value("driver");
    heat_fluids_driver_ =
      std::make_unique<SurrogateHeatDriver>(heat_comm.comm, heat_node);

    // Discover the ranks that are in each comm
    neutronics_ranks_ = gather_subcomm_ranks(comm_, neutronics_comm);
    heat_ranks_ = gather_subcomm_ranks(comm_, heat_comm);

    // Send rank of neutronics root to all procs
    neutronics_root_ = this->get_neutronics_driver().comm_.is_root() ? comm_.rank : -1;
    MPI_Allreduce(MPI_IN_PLACE, &neutronics_root_, 1, MPI_INT, MPI_MAX, comm_.comm);

    // Send rank of heat root to all procs
    heat_root_ = this->get_heat_driver().comm_.is_root() ? comm_.rank : -1;
    MPI_Allreduce(MPI_IN_PLACE, &heat_root_, 1, MPI_INT, MPI_MAX, comm_.comm);

    comm_report();

    init_tallies();
  }

  void comm_report()
  {
    char c[_POSIX_HOST_NAME_MAX];
    gethostname(c, _POSIX_HOST_NAME_MAX);
    std::string hostname{c};

    Comm world(MPI_COMM_WORLD);

    // Padding for fields
    int hostw = std::max(8UL, hostname.size()) + 2;
    int rankw = 7;

    for (int i = 0; i < world.size; ++i) {
      if (world.rank == i) {
        if (i == 0) {
          std::cout << "[ENRICO]: Communicator layout: " << std::endl;
          std::cout << std::left << std::setw(hostw) << "Hostname" << std::right
                    << std::setw(rankw) << "World" << std::right << std::setw(rankw)
                    << "Coup" << std::right << std::setw(rankw) << "Neut" << std::right
                    << std::setw(rankw) << "Heat" << std::endl;
        }
        std::cout << std::left << std::setw(hostw) << hostname << std::right
                  << std::setw(rankw) << world.rank << std::right << std::setw(rankw)
                  << comm_.rank << std::right << std::setw(rankw)
                  << this->get_neutronics_driver().comm_.rank << std::right
                  << std::setw(rankw) << this->get_heat_driver().comm_.rank << std::endl;
      }
      MPI_Barrier(world.comm);
    }
  }

  void init_tallies()
  {
    comm_.message("Initializing tallies");

    auto& neutronics = this->get_neutronics_driver();
    if (neutronics.active()) {
      neutronics.create_tallies();
    }
  }

  void execute()
  {
    auto& neutronics = get_neutronics_driver();
    MPI_Barrier(MPI_COMM_WORLD);
    comm_.message("Begin init_step");
    if (neutronics.active()) neutronics.init_step();
    MPI_Barrier(MPI_COMM_WORLD);
    comm_.message("Begin solve_step");
    if (neutronics.active()) neutronics.solve_step();
    MPI_Barrier(MPI_COMM_WORLD);
    comm_.message("Begin write_step");
    if (neutronics.active()) neutronics.write_step(i_timestep_, i_picard_);
    MPI_Barrier(MPI_COMM_WORLD);
    comm_.message("Begin finalize_step");
    if (neutronics.active()) neutronics.finalize_step();
  }

  NeutronicsDriver& get_neutronics_driver() const { return *neutronics_driver_; }
  HeatFluidsDriver& get_heat_driver() const { return *heat_fluids_driver_; }

  Comm comm_;
  int i_timestep_ = 0;
  int i_picard_ = 0;

  std::unique_ptr<NeutronicsDriver> neutronics_driver_;  //!< The neutronics driver
  std::unique_ptr<HeatFluidsDriver> heat_fluids_driver_; //!< The heat-fluids driver

  std::vector<int> heat_ranks_;
  std::vector<int> neutronics_ranks_;

  int neutronics_root_ = MPI_PROC_NULL;
  int heat_root_ = MPI_PROC_NULL;
};

}

int main(int argc, char* argv[])
{
  // Initialize MPI
  MPI_Init(&argc, &argv);
  enrico::init_mpi_datatypes();

  // Parse enrico.xml file
  pugi::xml_document doc;
  auto result = doc.load_file("enrico.xml");
  if (!result) {
    throw std::runtime_error{"Unable to load enrico.xml file"};
  }

  enrico::TestDriver driver{MPI_COMM_WORLD, doc.document_element()};

  driver.execute();

  enrico::free_mpi_datatypes();
  MPI_Finalize();
  return 0;
}
