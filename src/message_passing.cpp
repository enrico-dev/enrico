#include "enrico/message_passing.h"

#include "enrico/geom.h"

#include <mpi.h>

namespace enrico {

//==============================================================================
// Global variables
//==============================================================================

MPI_Datatype position_mpi_datatype{MPI_DATATYPE_NULL};

//==============================================================================
// Functions
//==============================================================================

void get_node_comms(MPI_Comm super_comm,
                    int procs_per_node,
                    MPI_Comm* sub_comm,
                    MPI_Comm* intranode_comm)
{

  // super_comm_rank is used as the "key" to retain ordering in the comm splits.
  // This can allow the sub_comm to retain some intent from the super_comm's proc layout
  int super_comm_rank;
  MPI_Comm_rank(super_comm, &super_comm_rank);

  // intranode_comm is an intermediate object.  It is only used to get an
  // intranode_comm_rank, which is used as the "color" in the final comm split.
  MPI_Comm_split_type(
    super_comm, MPI_COMM_TYPE_SHARED, super_comm_rank, MPI_INFO_NULL, intranode_comm);
  int intranode_comm_rank;
  MPI_Comm_rank(*intranode_comm, &intranode_comm_rank);

  // Finally, split the specified number of procs_per_node from the super_comm
  // We only want the comm where color == 0.  The second comm is destroyed.
  int color = intranode_comm_rank < procs_per_node ? 0 : 1;
  MPI_Comm_split(super_comm, color, super_comm_rank, sub_comm);
  if (color != 0)
    MPI_Comm_free(sub_comm);
}

void init_mpi_datatypes()
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
  MPI_Type_create_struct(3, blockcounts, displs, types, &position_mpi_datatype);
  MPI_Type_commit(&position_mpi_datatype);
}

void free_mpi_datatypes()
{
  MPI_Type_free(&position_mpi_datatype);
}

// Traits for mapping plain types to corresponding MPI types (ints)
template<>
MPI_Datatype get_mpi_type<char>()
{
  return MPI_CHAR;
}
template<>
MPI_Datatype get_mpi_type<short>()
{
  return MPI_SHORT;
}
template<>
MPI_Datatype get_mpi_type<int>()
{
  return MPI_INT;
}
template<>
MPI_Datatype get_mpi_type<long>()
{
  return MPI_LONG;
}
template<>
MPI_Datatype get_mpi_type<long long>()
{
  return MPI_LONG_LONG;
}
template<>
MPI_Datatype get_mpi_type<unsigned int>()
{
  return MPI_UNSIGNED;
}
template<>
MPI_Datatype get_mpi_type<unsigned long>()
{
  return MPI_UNSIGNED_LONG;
}
template<>
MPI_Datatype get_mpi_type<unsigned long long>()
{
  return MPI_UNSIGNED_LONG_LONG;
}

// Traits for mapping plain types to corresponding MPI types (reals)
template<>
MPI_Datatype get_mpi_type<float>()
{
  return MPI_FLOAT;
}
template<>
MPI_Datatype get_mpi_type<double>()
{
  return MPI_DOUBLE;
}

// Traits for mapping plain types to corresponding MPI types (user defined)
template<>
MPI_Datatype get_mpi_type<Position>()
{
  return position_mpi_datatype;
}

} // namespace enrico
