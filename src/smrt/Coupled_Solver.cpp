
#include <iostream>

#include "smrt/Coupled_Solver.h"
#include "stream/nek_driver.h"

namespace stream
{

// Constructor
Coupled_Solver::Coupled_Solver(std::shared_ptr<Assembly_Model> assembly,
                               const std::vector<double>&      z_edges,
                               const std::string&              shift_filename,
                               double                          power_norm,
                               MPI_Comm                        neutronics_comm,
                               MPI_Comm                        th_comm)
    : d_power_norm(power_norm)
{
    d_shift_solver = std::make_shared<stream::Shift_Solver>(
        assembly,
        shift_filename,
        z_edges);

    d_nek_solver = std::make_shared<NekDriver>(th_comm);

    d_th_num_local  = d_nek_solver->nelt_;
    d_th_num_global = d_nek_solver->nelgt_;

    this->init_mpi_datatypes();

    // Allocate fields (on global T/H mesh for now)
    d_temperatures.resize(d_th_num_global, 565.0);
    d_densities.resize(d_th_num_global, 0.75);
    d_powers.resize(d_th_num_global, 0.0);

    std::vector<Position> local_centroids(d_th_num_local);
    std::vector<double> local_volumes(d_th_num_local);
    for (int elem = 0; elem < d_th_num_local; ++elem)
    {
        local_centroids[elem] = d_nek_solver->get_local_elem_centroid(elem+1);
        local_volumes[elem]   = d_nek_solver->get_local_elem_volume(elem+1);
        assert(!std::isnan(local_centroids[elem].x));
        assert(!std::isnan(local_centroids[elem].y));
        assert(!std::isnan(local_centroids[elem].z));
        assert(local_volumes[elem] > 0.0);
    }
    auto global_centroids = this->local_to_global(local_centroids);
    auto global_volumes   = this->local_to_global(local_volumes);

    // Sanity check on centroids and volumes
    for (const auto& c : global_centroids)
    {
        assert(!std::isnan(c.x));
        assert(!std::isnan(c.y));
        assert(!std::isnan(c.z));
    }
    for (const auto& v : global_volumes)
        assert(v > 0.0);

    // Register centroids and volumes with Shift
    d_shift_solver->set_centroids_and_volumes(
        global_centroids,
        global_volumes);
}

// Destructor
Coupled_Solver::~Coupled_Solver()
{
}

// Solve coupled problem by iterating between neutronics and T/H
void Coupled_Solver::solve()
{
    d_nek_solver->solve_step();

    d_shift_solver->solve(d_temperatures, d_densities, d_powers);

    this->free_mpi_datatypes();
}

//
// Private Implementation
//

// Set up MPI datatype
// Currently, this sets up only position_mpi_datatype
void Coupled_Solver::init_mpi_datatypes()
{
    Position p;
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
    MPI_Type_create_struct(3, blockcounts, displs, types, &d_position_mpi_type);
    MPI_Type_commit(&d_position_mpi_type);
}

// Free user-defined MPI types
void Coupled_Solver::free_mpi_datatypes()
{
    MPI_Type_free(&d_position_mpi_type);
}

// Traits for mapping plain types to corresponding MPI types
template <>
MPI_Datatype Coupled_Solver::get_mpi_type<double>() const
{
    return MPI_DOUBLE;
}
template <>
MPI_Datatype Coupled_Solver::get_mpi_type<Position>() const
{
    return d_position_mpi_type;
}

// Gather local distributed field into global replicated field
template <typename T>
std::vector<T> Coupled_Solver::local_to_global(
    const std::vector<T>& local_field) const
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

    th_comm.Bcast(global_field.data(),
                  d_th_num_global,
                  get_mpi_type<T>());

    return global_field;
}

// Scatter global replicated field into local distributed field
template <typename T>
std::vector<T> Coupled_Solver::global_to_local(
    const std::vector<T>& global_field) const
{
    const auto& th_comm = d_nek_solver->comm_;

    // Local elements are contiguous in global array with rank-specific offset
    int offset = d_nek_solver->local_displs_[offset];
    std::vector<T> local_field(d_th_num_local);
    for (int elem = 0; elem < d_th_num_local; ++elem)
        local_field[elem] = global_field[elem + offset];

    return local_field;
}

} // end namespace stream

