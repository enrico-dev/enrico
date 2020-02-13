//! \file message_passing.h
//! Utility functions for constucting MPI communicators
#ifndef ENRICO_MPI_TYPES_H
#define ENRICO_MPI_TYPES_H

#include <mpi.h>

//! The ENRICO namespace
namespace enrico {

//==============================================================================
// Global variables
//==============================================================================

extern MPI_Datatype position_mpi_datatype;

//==============================================================================
// Functions
//==============================================================================

//! Create MPI datatype for Position struct
void init_mpi_datatypes();

//! Free any MPI datatypes
void free_mpi_datatypes();

//! Map types to corresponding MPI datatypes
template<typename T>
MPI_Datatype get_mpi_type();

} // namespace enrico

#endif // ENRICO_MPI_TYPES_H
