#include "enrico/mpi_types.h"

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

template<>
MPI_Datatype get_mpi_type<bool>()
{
  return MPI_CXX_BOOL;
}

// Traits for mapping plain types to corresponding MPI types (user defined)
template<>
MPI_Datatype get_mpi_type<Position>()
{
  return position_mpi_datatype;
}

} // namespace enrico
