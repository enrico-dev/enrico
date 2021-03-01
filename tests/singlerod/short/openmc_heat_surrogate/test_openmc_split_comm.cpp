#include "enrico/comm_split.h"
#include "enrico/openmc_driver.h"
#include "enrico/surrogate_heat_driver.h"
#include <iomanip>
#include <mpi.h>
#include <unistd.h>

namespace enrico {

class TestDriver {
public:
  enum class Initial { neutronics, heat };

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

    init_mappings();
    init_tallies();
    init_volumes();

    check_volumes();

    init_fluid_mask();

    init_temperatures();
    init_densities();
    init_heat_source();
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

  void init_mappings()
  {
    comm_.message("Initializing mappings");
    const auto& heat = this->get_heat_driver();
    auto& neutronics = this->get_neutronics_driver();

    // Send and recv buffers
    std::vector<Position> centroids_send, centroids_recv;
    decltype(elem_to_cell_) elem_to_cell_send;

    for (const auto& heat_rank : heat_ranks_) {
      // Set mapping of local elem --> global cell handle
      if (comm_.rank == heat_rank) {
        centroids_send = heat.centroid_local();
      }
      this->comm_.send_and_recv(
        centroids_recv, neutronics_root_, centroids_send, heat_rank);
      if (comm_.rank == neutronics_root_) {
        elem_to_cell_send = neutronics.find(centroids_recv);
      }
      this->comm_.send_and_recv(
        elem_to_cell_, heat_rank, elem_to_cell_send, neutronics_root_);
      comm_.Barrier();
    }
    if (heat.active()) {
      // Set mapping of global cell handle -> local elem
      for (gsl::index e = 0; e < elem_to_cell_.size(); ++e) {
        auto c = elem_to_cell_[e];
        cell_to_elems_[c].push_back(e); // Use [ ] instead of at() to insert new item
      }
      // Get list of all global cell handles
      for (const auto& kv : cell_to_elems_) {
        cells_.push_back(kv.first);
      }
    }
  }

  void init_volumes()
  {
    comm_.message("Initializing volumes");
    const auto& heat = this->get_heat_driver();
    const auto& neutronics = this->get_neutronics_driver();

    if (heat.active()) {
      elem_volumes_ = heat.volume_local();
      for (const auto& c : cells_) {
        double V = 0.0;
        for (const auto& e : cell_to_elems_.at(c)) {
          V += elem_volumes_.at(e);
        }
        cell_volumes_.push_back(V);
      }
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

  void check_volumes()
  {
    comm_.message("Initializing volumes");
    const auto& neutronics = this->get_neutronics_driver();

    // An array of global cell volumes, which will be accumulated from local cell volumes.
    std::map<CellHandle, double> glob_volumes;

    // Get all local cell volumes from heat ranks and sum them into the global cell
    // volumes.
    for (const auto& heat_rank : heat_ranks_) {
      decltype(cells_) cells_recv;
      decltype(cell_volumes_) cell_volumes_recv;
      comm_.send_and_recv(cells_recv, neutronics_root_, cells_, heat_rank);
      comm_.send_and_recv(cell_volumes_recv, neutronics_root_, cell_volumes_, heat_rank);
      if (comm_.rank == neutronics_root_) {
        for (gsl::index i = 0; i < cells_recv.size(); ++i) {
          glob_volumes[cells_recv.at(i)] += cell_volumes_recv.at(i);
        }
      }
    }
    comm_.Barrier();

    if (comm_.rank == neutronics_root_) {
      // Compare volume from neutron driver to accumulated volume
      for (const auto& kv : glob_volumes) {
        auto v_neutronics = neutronics.get_volume(kv.first);
        std::stringstream msg;
        msg << "Cell " << neutronics.cell_label(kv.first) << ", V = " << v_neutronics
            << " (Neutronics), " << kv.second << " (Accumulated from Heat/Fluids)";
        comm_.message(msg.str());
      }
    }
    comm_.Barrier();
  }

  void init_fluid_mask()
  {
    comm_.message("Initializing cell fluid mask");
    auto& heat = this->get_heat_driver();

    if (heat.active()) {
      auto elem_fluid_mask = heat.fluid_mask_local();
      for (const auto& kv : cell_to_elems_) {

        auto& elems = kv.second;
        auto in_fluid = elem_fluid_mask.at(elems.at(0));
        for (gsl::index i = 1; i < elems.size(); ++i) {
          if (in_fluid != elem_fluid_mask.at(elems.at(i))) {
            throw std::runtime_error("ENRICO detected a neutronics cell contains both "
                                     "fluid and solid T/H elements.");
          }
        }
        cell_fluid_mask_.push_back(in_fluid);
      }
    }
  }

  void init_temperatures()
  {
    comm_.message("Initializing temperatures");
    const auto& neutronics = this->get_neutronics_driver();
    const auto& heat = this->get_heat_driver();

    // Every heat rank keeps track of its own local cell T
    if (heat.active()) {
      auto sz = static_cast<unsigned long>(cells_.size());
      cell_temperatures_.resize({sz});
      cell_temperatures_prev_.resize({sz});
    }

    if (temperature_ic_ == Initial::neutronics) {
      // Send and recv buffers
      decltype(cells_) cells_recv;
      decltype(cell_temperatures_) cell_temperatures_send;
      // The neutronics root sends cell T to each heat rank
      for (const auto& heat_rank : heat_ranks_) {
        comm_.send_and_recv(cells_recv, neutronics_root_, cells_, heat_rank);
        if (comm_.rank == neutronics_root_) {
          const auto sz = static_cast<unsigned long>(cells_recv.size());
          cell_temperatures_send.resize({sz});
          for (gsl::index i = 0; i < sz; ++i) {
            cell_temperatures_send.at(i) = neutronics.get_temperature(cells_recv.at(i));
          }
        }
        comm_.send_and_recv(
          cell_temperatures_, heat_rank, cell_temperatures_send, neutronics_root_);
      }
    } else if (temperature_ic_ == Initial::heat) {
      // * This sets l_cel_temps_ on the the coupling_root, based on the
      //   temperatures received from the heat solver.
      // * We do not want to apply underrelaxation here (and at this point,
      //   there is no previous iterate of temperature, anyway).
      update_temperature(false);
    }

    // In both cases, only temperatures_ was set, so we explicitly set temperatures_prev_
    if (heat.active()) {
      std::copy(cell_temperatures_.begin(),
                cell_temperatures_.end(),
                cell_temperatures_prev_.begin());
    }
  }

  void update_temperature(bool relax)
  {
    comm_.message("Updating temperature");
    auto& neutronics = this->get_neutronics_driver();
    auto& heat = this->get_heat_driver();

    // Step 1: On heat ranks, assign the current iterate of T to the previous
    if (relax && heat.active()) {
      std::copy(cell_temperatures_.begin(),
                cell_temperatures_.end(),
                cell_temperatures_prev_.begin());
    }

    // Step 2: On heat ranks, compute cell T
    if (heat.active()) {
      auto elem_temperatures = heat.temperature_local();
      for (gsl::index i = 0; i < cells_.size(); ++i) {
        double T_avg = 0.0;
        double V_tot = 0.0;
        for (const auto& e : cell_to_elems_.at(cells_.at(i))) {
          double T = elem_temperatures.at(e);
          double V = elem_volumes_.at(e);
          T_avg += T * V;
          V_tot += V;
        }
        T_avg /= V_tot;
        Ensures(T_avg > 0.0);
        cell_temperatures_.at(i) = T_avg;
      }
      // Apply relaxation to local cell temperatures
      if (relax) {
        if (alpha_T_ == ROBBINS_MONRO) {
          int n = i_picard_ + 1;
          cell_temperatures_ =
            cell_temperatures_ / n + (1. - 1. / n) * cell_temperatures_prev_;
        } else {
          cell_temperatures_ =
            alpha_T_ * cell_temperatures_ + (1.0 - alpha_T_) * cell_temperatures_prev_;
        }
      }
    }

    // Step 3: On neutron ranks, accumulate cell volumes from all heat ranks
    std::map<CellHandle, double> T_dot_V;
    std::map<CellHandle, double> cell_V;
    decltype(cells_) cells_recv;
    decltype(cell_volumes_) cell_volumes_recv;
    decltype(cell_temperatures_) cell_temperatures_recv;
    for (const auto& heat_rank : heat_ranks_) {
      comm_.send_and_recv(cells_recv, neutronics_root_, cells_, heat_rank);
      neutronics.comm_.broadcast(cells_recv);

      comm_.send_and_recv(
        cell_temperatures_recv, neutronics_root_, cell_temperatures_, heat_rank);
      neutronics.comm_.broadcast(cell_temperatures_recv);

      comm_.send_and_recv(cell_volumes_recv, neutronics_root_, cell_volumes_, heat_rank);
      neutronics.comm_.broadcast(cell_volumes_recv);

      if (neutronics.active()) {
        for (gsl::index i = 0; i < cells_recv.size(); ++i) {
          auto T = cell_temperatures_recv.at(i);
          auto V = cell_volumes_recv.at(i);
          cell_V[cells_recv.at(i)] += V;
          T_dot_V[cells_recv.at(i)] += T * V;
        }
      }
    }
    for (const auto& kv : T_dot_V) {
      neutronics.set_temperature(kv.first, kv.second / cell_V.at(kv.first));
    }
  }

  void update_density(bool relax)
  {
    comm_.message("Updating density");
    auto& neutronics = this->get_neutronics_driver();
    auto& heat = this->get_heat_driver();

    // Step 1: On heat ranks, assign the current iterate of rho to the previous
    if (relax && heat.active()) {
      std::copy(
        cell_densities_.cbegin(), cell_densities_.cend(), cell_densities_.begin());
    }

    // Step 2: On heat ranks, get cell rho
    if (heat.active()) {
      auto elem_densities = heat.density_local();

      for (gsl::index i = 0; i < cells_.size(); ++i) {
        if (cell_fluid_mask_.at(i) == 1) {
          double rho_avg = 0.0;
          double V_tot = 0.0;
          for (const auto& e : cell_to_elems_.at(cells_.at(i))) {
            rho_avg += elem_densities.at(e) * elem_volumes_.at(e);
            V_tot += elem_volumes_.at(e);
          }
          rho_avg /= V_tot;
          Ensures(rho_avg > 0.0);
          cell_densities_.at(i) = rho_avg;
        }
      }
      if (relax) {
        if (alpha_rho_ == ROBBINS_MONRO) {
          int n = i_picard_ + 1;
          cell_densities_ = cell_densities_ / n + (1. - 1. / n) * cell_densities_prev_;
        } else {
          cell_densities_ =
            alpha_rho_ * cell_densities_ + (1.0 - alpha_rho_) * cell_densities_prev_;
        }
      }
    }

    // Step 3: On neutron ranks, accumulate cell rho from all heat ranks
    std::map<CellHandle, double> rho_dot_V;
    std::map<CellHandle, double> cell_V;
    decltype(cells_) cells_recv;
    decltype(cell_volumes_) cell_volumes_recv;
    decltype(cell_densities_) cell_densities_recv;
    decltype(cell_fluid_mask_) cell_fluid_mask_recv;

    for (const auto& heat_rank : heat_ranks_) {
      comm_.send_and_recv(cells_recv, neutronics_root_, cells_, heat_rank);
      neutronics.comm_.broadcast(cells_recv);

      comm_.send_and_recv(cell_volumes_recv, neutronics_root_, cell_volumes_, heat_rank);
      neutronics.comm_.broadcast(cell_volumes_recv);

      comm_.send_and_recv(
        cell_densities_recv, neutronics_root_, cell_densities_, heat_rank);
      neutronics.comm_.broadcast(cell_densities_recv);

      comm_.send_and_recv(
        cell_fluid_mask_recv, neutronics_root_, cell_fluid_mask_, heat_rank);

      if (neutronics.active()) {
        for (gsl::index i = 0; i < cells_recv.size(); ++i) {
          if (cell_fluid_mask_recv.at(i) == 1) {
            auto rho = cell_densities_recv.at(i);
            auto V = cell_volumes_recv.at(i);
            cell_V[cells_recv.at(i)] += V;
            rho_dot_V[cells_recv.at(i)] += rho * V;
          }
        }
      }
    }
    for (const auto& kv : rho_dot_V) {
      neutronics.set_density(kv.first, kv.second / cell_V.at(kv.first));
    }
  }

  void init_densities()
  {
    comm_.message("Initializing densities");
    const auto& neutronics = this->get_neutronics_driver();
    const auto& heat = this->get_heat_driver();

    if (heat.active()) {
      auto sz = static_cast<unsigned long>(cells_.size());
      cell_densities_.resize({sz});
      cell_densities_prev_.resize({sz});
    }

    if (density_ic_ == Initial::neutronics) {
      decltype(cells_) cells_recv;
      decltype(cell_densities_) cell_densities_send;
      for (const auto& heat_rank : heat_ranks_) {
        comm_.send_and_recv(cells_recv, neutronics_root_, cells_, heat_rank);
        if (comm_.rank == neutronics_root_) {
          const auto sz = static_cast<unsigned long>(cells_recv.size());
          cell_densities_send.resize({sz});
          for (gsl::index i = 0; i < sz; ++i) {
            cell_densities_send.at(i) = neutronics.get_density(cells_recv.at(i));
          }
        }
        comm_.send_and_recv(
          cell_densities_, heat_rank, cell_densities_send, neutronics_root_);
      }
    } else if (density_ic_ == Initial::heat) {
      // * We do not want to apply underrelaxation here (and at this point,
      //   there is no previous iterate of density, anyway).
      update_density(false);
    }

    // In both cases, we need to explicitly set densities_prev_
    if (heat.active()) {
      std::copy(
        cell_densities_.cbegin(), cell_densities_.cend(), cell_densities_prev_.begin());
    }
  }

  void init_heat_source()
  {
    comm_.message("Initializing heat source");

    if (this->heat_fluids_driver_->active()) {
      auto sz = {cells_.size()};
      cell_heat_ = xt::empty<double>(sz);
      cell_heat_prev_ = xt::empty<double>(sz);
    }
  }

  void execute()
  {
    auto& neutronics = get_neutronics_driver();
    neutronics.init_step();
    neutronics.solve_step();
    // neutronics.write_step(i_timestep_, i_picard_);
    neutronics.finalize_step();
  }

  NeutronicsDriver& get_neutronics_driver() const { return *neutronics_driver_; }
  HeatFluidsDriver& get_heat_driver() const { return *heat_fluids_driver_; }

  constexpr static double ROBBINS_MONRO = -1.0;

  Comm comm_;
  int i_timestep_ = 0;
  int i_picard_ = 0;

  double alpha_{1.0};
  double alpha_T_{alpha_};
  double alpha_rho_{alpha_};

  Initial temperature_ic_{Initial::neutronics};
  Initial density_ic_{Initial::neutronics};

  std::unique_ptr<NeutronicsDriver> neutronics_driver_;  //!< The neutronics driver
  std::unique_ptr<HeatFluidsDriver> heat_fluids_driver_; //!< The heat-fluids driver

  std::vector<int> heat_ranks_;
  std::vector<int> neutronics_ranks_;

  int neutronics_root_ = MPI_PROC_NULL;
  int heat_root_ = MPI_PROC_NULL;

  xt::xtensor<double, 1> cell_densities_;
  xt::xtensor<double, 1> cell_densities_prev_;
  xt::xtensor<double, 1> cell_temperatures_;
  xt::xtensor<double, 1> cell_temperatures_prev_;
  xt::xtensor<double, 1> cell_heat_;
  xt::xtensor<double, 1> cell_heat_prev_;
  std::vector<int> cell_fluid_mask_;
  std::vector<CellHandle> elem_to_cell_;
  std::vector<CellHandle> cells_;
  std::map<CellHandle, std::vector<int32_t>> cell_to_elems_;
  std::vector<double> cell_volumes_;
  std::vector<double> elem_volumes_;
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
