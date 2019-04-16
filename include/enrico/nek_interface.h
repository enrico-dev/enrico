//! \file nek_interface.h
//! Functions for accessing Nek5000 routines and data structures
#ifndef ENRICO_NEK_INTERFACE_H
#define ENRICO_NEK_INTERFACE_H

#include "geom.h"
#include "nek_mangling.h"

// TODO: The Doxygen doesn't output these global functions

extern "C" {

//!  One-time initialization of Nek5000
//!
//!  This is the name-mangled version of nek_init from libnek5000
//!
//!  \param intracomm An exisiting MPI communicator used to run Nek5000
void C2F_nek_init(const int* intracomm);

//! One-time finalization of Nek5000
//!
//! This is the name-mangled version of nek_end from libnek5000
void C2F_nek_end();

//! Do an entire solve (all timesteps) of Nek5000
//!
//! This is the name-mangled version of nek_solve from libnek5000
void C2F_nek_solve();

//! Reset the counters necessary to resume timestepping at the next Picard iteration
void nek_reset_counters();

//! Get the coordinates of a global element's centroid
//!
//! The returned coordinate is dimensionless.  Its units depend on the unit system that
//! was used to setup the Nek5000 problem. The user must handle any necessary
//! conversions.
//!
//! \param[in] global_elem A global element ID
//! \param[out] centroid The **dimensionless* position of the global element's centroid
//! \return Error value
int nek_get_global_elem_centroid(int global_elem, enrico::Position* centroid);

//! Get the coordinates of a local element's centroid
//!
//! The returned coordinate is dimensionless.  Its units depend on the unit system that
//! was used to setup the Nek5000 problem. The user must handle any necessary
//! conversions.
//!
//! \param[in] local_elem A local element ID
//! \param[out] centroid The **dimensionless* position of the local element's centroid
//! \return Error value
int nek_get_local_elem_centroid(int local_elem, enrico::Position* centroid);

//! Get the volume of a local element
//!
//! The returned volume is dimensionless.  Its units depend on the unit system that was
//! used to setup the Nek5000 problem. The user must handle any necessary conversions.
//!
//! \param local_elem  A local element ID
//! \param volume  The **dimensionless** volume of the local element
//! \return  Error value
int nek_get_local_elem_volume(int local_elem, double* volume);

//! Get the volume-averaged temperature of a local element
//!
//! The returned temperature is dimensionless.  Its units depend on the unit system that
//! was used to setup the Nek5000 problem. The user must handle any necessary
//! conversions.
//!
//! \param local_elem  A local element ID
//! \param temperature  The **dimensionless** volume-averaged temperature of the local
//! element \return  Error value
int nek_get_local_elem_temperature(int local_elem, double* temperature);

//! Get the global element ID for a given local element
//! \param local_elem  A local element ID
//! \return The corresponding global element ID
int nek_get_global_elem(int local_elem);

//! Get the local element ID for a given global element
//! \param global_elem  A global element ID
//! \return The corresponding local element ID
int nek_get_local_elem(int global_elem);

//! Get lelg, the maximum number of global elements
//!
//! For Nek5000, lelg is the compile-time constant for the maximum number of global
//! elements.  At runtime, the user will choose an actual number of global elements
//! (nelg) for their problem, which must be <= lelg.
//!
//! \return  The upper bound on number of global elements
int nek_get_lelg();

//! Get lelt, the maximum number of local elements
//!
//! For Nek5000, lelt is the compile-time constant for the maximum number of local
//! elements.  At runtime, the user will choose an actual number of local elements
//! (nelg) for their problem, which must be <= lelt.
//!
//! \return  The upper bound on number of local elements
int nek_get_lelt();

//! Get lx1, the maximum number of GLL gridpoints in the x-direction
//!
//! For Nek5000, lx1 is the compile-time constant for the maximum number GLL gridpoints
//! in the x-direction.  At runtime, the user will choose an actual number of GLL
//! gridpoints for their problem, which must be <= lx1.
//!
//! \return  The upper bound on number of GLL gridpoints in the x-direction
int nek_get_lx1();

//! Get lelt, the number of local elements
//!
//! \return  The number of local elements
int nek_get_nelt();

//! Get nelgt, the number of global elements
//! \return  nelgt The number of global elements
int nek_get_nelgt();

//! Return true if a global element is in a given MPI rank
//! \param A global element ID
//! \param An MPI rank
//! \return True if the global element ID is in the given rank
int nek_global_elem_is_in_rank(int global_elem, int rank);

//! Return true if a local element is in the fluid region
//! \param local_elem  A local element ID
//! \return 1 if the local element is in fluid; 0 otherwise
int nek_local_elem_is_in_fluid(int local_elem);

//! Return true if a global element is in the fluid region
//! \param global_elem  A global element ID
//! \return 1 if the global element is in fluid; 0 otherwise
int nek_global_elem_is_in_fluid(int global_elem);

//! Set the heat source for a given local element
//!
//! The units of heat must match on the unit system that was used to setup the Nek5000
//! problem (presumably W/cm^3). The caller must handle any necessary conversions.
//!
//! \param local_elem A local element ID
//! \param heat A heat source term
//! \return Error code
int nek_set_heat_source(int local_elem, double heat);

//! Unimplemented.  C2F_nek_init is used instead.
void nek_init_step();

//! Unimplemented. C2F_nek_solve is used instead.
void nek_step();

//! Unimplemented. C2F_nek_end is used instead.
void nek_finalize_step();
};

#endif // ENRICO_NEK_INTERFACE_H
