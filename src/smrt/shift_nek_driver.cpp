#include <iostream>

#include "enrico/error.h"
#include "enrico/nek_driver.h"
#include "nek5000/core/nek_interface.h"
#include "smrt/shift_nek_driver.h"

namespace enrico {

// Constructor
ShiftNekDriver::ShiftNekDriver(std::shared_ptr<Assembly_Model> assembly,
                               const std::vector<double>& z_edges,
                               const std::string& shift_filename,
                               MPI_Comm neutronics_comm,
                               MPI_Comm th_comm) :
  SmrtCoupledDriver()
{
  d_shift_solver =
    std::make_shared<enrico::ShiftDriver>(assembly, shift_filename, z_edges);

  // Build Nek driver
  {
    // TODO: Belongs to main.cpp
    pugi::xml_document doc;
    auto result = doc.load_file("enrico.xml");
    if (!result) {
      throw std::runtime_error{"Unable to load enrico.xml file"};
    }

    // Get root element
    auto root = doc.document_element();

    // TODO: Belongs to CoupledDriver base class
    power_ = root.child("power").text().as_double();
    max_picard_iter_ = root.child("max_picard_iter").text().as_int();

    double pressure_bc = root.child("pressure_bc").text().as_double();

    d_nek_solver = std::make_shared<NekDriver>(
        th_comm, pressure_bc, root.child("nek5000"));
  }

  d_th_num_local = d_nek_solver->nelt_;
  d_th_num_global = d_nek_solver->nelgt_;

  this->init_mpi_datatypes();

  // Allocate fields (on global T/H mesh for now)
  d_temperatures.resize(d_th_num_local, 565.0);
  d_densities.resize(d_th_num_local, 0.75);
  d_powers.resize(d_th_num_local, power_);

  std::vector<Position> local_centroids(d_th_num_local);
  std::vector<double> local_volumes(d_th_num_local);
  for (int elem = 0; elem < d_th_num_local; ++elem) {
    local_centroids[elem] = d_nek_solver->centroid_at(elem + 1);
    local_volumes[elem] = d_nek_solver->volume_at(elem + 1);
    assert(!std::isnan(local_centroids[elem].x));
    assert(!std::isnan(local_centroids[elem].y));
    assert(!std::isnan(local_centroids[elem].z));
    assert(local_volumes[elem] > 0.0);
  }
  auto global_centroids = this->local_to_global(local_centroids);
  auto global_volumes = this->local_to_global(local_volumes);

  // Sanity check on centroids and volumes
  for (const auto& c : global_centroids) {
    assert(!std::isnan(c.x));
    assert(!std::isnan(c.y));
    assert(!std::isnan(c.z));
  }
  for (const auto& v : global_volumes)
    assert(v > 0.0);

  // Register centroids and volumes with Shift
  d_shift_solver->set_centroids_and_volumes(local_centroids, local_volumes);

  for (int rank = 0; rank < nemesis::nodes(); ++rank) {
    if (rank == nemesis::node()) {
      std::cout << "Element counts on " << rank << ": ";
      for (auto val : d_nek_solver->local_counts_)
        std::cout << val << " ";
      std::cout << std::endl;
      std::cout << "Displacements on " << rank << ": ";
      for (auto val : d_nek_solver->local_displs_)
        std::cout << val << " ";
      std::cout << std::endl;
    }
    nemesis::global_barrier();
  }
}

// Destructor
ShiftNekDriver::~ShiftNekDriver() {}

// Solve coupled problem by iterating between neutronics and T/H
void ShiftNekDriver::solve()
{
  // Loop to convergence or fixed iteration count
  for (int iteration = 0; iteration < max_picard_iter_; ++iteration) {
    // Set heat source in Nek
    for (int elem = 0; elem < d_th_num_local; ++elem) {
      err_chk(d_nek_solver->set_heat_source_at(elem + 1, d_powers[elem]),
              "Error setting heat source for local element " + std::to_string(elem + 1));
    }

    // Solve Nek problem
    d_nek_solver->solve_step();

    // Get temperatures from Nek
    for (int elem = 0; elem < d_th_num_local; ++elem) {
      // Normalization for incorrect Gauss point averaging
      constexpr double nek_normalization = 1.0 / 200.0;
      d_temperatures[elem] =
        d_nek_solver->temperature_at(elem + 1) * nek_normalization;
    }

    for (int rank = 0; rank < nemesis::nodes(); ++rank) {
      if (rank == nemesis::node()) {
        std::cout << "Temperature on " << rank << ": ";
        for (auto val : d_temperatures)
          std::cout << val << " ";
        std::cout << std::endl;
      }
      nemesis::global_barrier();
    }

    // Solve Shift problem
    d_shift_solver->solve(d_temperatures, d_densities, d_powers);

    // Apply power normalization
    normalize_power();

    for (int rank = 0; rank < nemesis::nodes(); ++rank) {
      if (rank == nemesis::node()) {
        std::cout << "Power on " << rank << ": ";
        for (auto val : d_powers)
          std::cout << val << " ";
        std::cout << std::endl;
      }
      nemesis::global_barrier();
    }
  }

  this->free_mpi_datatypes();
}

//
// Private Implementation
//

// Apply power normalization
void ShiftNekDriver::normalize_power()
{
  double total_power = 0.0;
  for (int elem = 0; elem < d_th_num_local; ++elem) {
    total_power += d_powers[elem] * d_nek_solver->volume_at(elem + 1);
  }
  nemesis::global_sum(total_power);

  // Apply normalization factor
  double norm_factor = power_ / total_power;
  for (auto& val : d_powers)
    val *= norm_factor;
}

// Set up MPI datatype
// Currently, this sets up only position_mpi_datatype
void ShiftNekDriver::init_mpi_datatypes()
{
  d_position_mpi_type = define_position_mpi_type();
}

// Free user-defined MPI types
void ShiftNekDriver::free_mpi_datatypes()
{
  MPI_Type_free(&d_position_mpi_type);
}

// Traits for mapping plain types to corresponding MPI types
template<>
MPI_Datatype ShiftNekDriver::get_mpi_type<double>() const
{
  return MPI_DOUBLE;
}
template<>
MPI_Datatype ShiftNekDriver::get_mpi_type<Position>() const
{
  return d_position_mpi_type;
}

// Gather local distributed field into global replicated field
template<typename T>
std::vector<T> ShiftNekDriver::local_to_global(const std::vector<T>& local_field) const
{
  assert(local_field.size() == d_th_num_local);
  const auto& th_comm = d_nek_solver->comm_;
  std::vector<T> global_field(d_th_num_global);

  // Do a gather only root and then broadcast
  th_comm.Gatherv(local_field.data(),
                  d_th_num_local,
                  get_mpi_type<T>(),
                  global_field.data(),
                  d_nek_solver->local_counts_.data(),
                  d_nek_solver->local_displs_.data(),
                  get_mpi_type<T>());

  th_comm.Bcast(global_field.data(), d_th_num_global, get_mpi_type<T>());

  return global_field;
}

// Scatter global replicated field into local distributed field
template<typename T>
std::vector<T> ShiftNekDriver::global_to_local(const std::vector<T>& global_field) const
{
  const auto& th_comm = d_nek_solver->comm_;

  // Local elements are contiguous in global array with rank-specific offset
  int offset = d_nek_solver->local_displs_[offset];
  std::vector<T> local_field(d_th_num_local);
  for (int elem = 0; elem < d_th_num_local; ++elem)
    local_field[elem] = global_field[elem + offset];

  return local_field;
}

} // end namespace enrico
