//! \file nek_driver.h
//! Driver to initialze and run Nek5000 in stages
#ifndef ENRICO_NEK_DRIVER_H
#define ENRICO_NEK_DRIVER_H

#include "driver.h"
#include "geom.h"
#include "mpi.h"

#include <pugixml.hpp>
#include <string>
#include <vector>

namespace enrico {

//! Driver to initialze and run Nek5000 in stages.
class NekDriver : public Driver {
public:
  //! Initializes Nek5000 with the given MPI communicator.
  //!
  //! A wrapper for the nek_init() routine in Nek5000.  Also retreives common-block
  //! variables from Nek5000 and uses them to set the corresponding attributes in
  //! NekDriver.
  //!
  //! \param comm  The MPI communicator used to initialze Nek5000
  explicit NekDriver(MPI_Comm comm, pugi::xml_node xml_root);

  //! Finalizes Nek5000.
  //!
  //! A wrapper for the nek_end() routine in Nek5000.
  ~NekDriver();

  void init_session_name();

  //! Runs all timesteps for a heat/fluid solve in Nek5000.
  //!
  //! A wraper for the nek_solve() routine in libnek5000.  This includes the necessary
  //! initialization and finalization for each step.
  void solve_step() final;

  //! Get the coordinate of a global element's centroid.
  //!
  //! The coordinate is dimensionless.  Its units depend on the unit system used that was
  //! used to setup the Nek problem. The user must handle any necessary conversions.
  //!
  //! \param global_elem The global index of the desired element
  //! \return The dimensionless coordinate of the element's centroid
  Position get_global_elem_centroid(int global_elem) const;

  //! Get the coordinate of a local element's centroid.
  //!
  //! The coordinate is dimensionless.  Its units depend on the unit system used that was
  //! used to setup the Nek problem. The user must handle any necessary conversions.
  //!
  //! \param local_elem The local index of the desired element
  //! \return The dimensionless coordinate of the element's centroid
  Position get_local_elem_centroid(int local_elem) const;

  //! Get the volume of a local element
  //!
  //! The volume is dimensionless.  Its units depend on the unit system used that was used
  //! to setup the Nek problem. The user must handle any necessary conversions.
  //!
  //! \param local_elem The local index of the desired element
  //! \return The dimensionless Volume of the element
  double get_local_elem_volume(int local_elem) const;

  //! Get the volume-averaged temperature of a local element
  //!
  //! The returned temperature is dimensionless.  Its units depend on the unit system that
  //! was used to setup the Nek5000 problem. The user must handle any necessary
  //! conversions.
  //!
  //! \param local_elem A local element ID
  //! \return The volume-averaged temperature of the element
  double get_local_elem_temperature(int local_elem) const;

  //! Return true if a global element is in a given MPI rank
  //! \param A global element ID
  //! \param An MPI rank
  //! \return True if the global element ID is in the given rank
  bool global_elem_is_in_rank(int global_elem) const;

  //! Return true if a global element is in the fluid region
  //! \param global_elem  A global element ID
  //! \return 1 if the global element is in fluid; 0 otherwise
  int global_elem_is_in_fluid(int global_elem) const;

  //! Return true if a local element is in the fluid region
  //! \param local_elem  A local element ID
  //! \return 1 if the local element is in fluid; 0 otherwise
  int local_elem_is_in_fluid(int local_elem) const;

  //! Set the heat source for a given local element
  //!
  //! The units of heat must match on the unit system that was used to setup the Nek5000
  //! problem (presumably W/cm^3). The caller must handle any necessary conversions.
  //!
  //! \param local_elem A local element ID
  //! \param heat A heat source term
  //! \return Error code
  int set_heat_source(int local_elem, double heat) const;

  //! Initialize the counts and displacements of local elements for each MPI Rank.
  void init_displs();

  std::string casename_; //!< Nek5000 casename (name of .rea file)
  int lelg_;             //!< upper bound on number of mesh elements
  int lelt_;             //!< upper bound on number of mesh elements per rank
  int lx1_;              //!< polynomial order of the solution
  int nelgt_;            //!< total number of mesh elements
  int nelt_;             //!< number of local mesh elements

  //! The number of local elements in each rank.
  std::vector<int> local_displs_;

  //! The displacements of local elements, relative to rank 0. Used in an MPI gatherv
  //! operation.
  std::vector<int> local_counts_;

  // Intended to be the local-to-global element ordering, as ensured by a Gatherv
  // operation. It is currently unused, as the coupling does not need to know the
  // local-global ordering.
  // std::vector<int> local_ordering_;
};

} // namespace enrico

#endif // ENRICO_NEK_DRIVER_H
