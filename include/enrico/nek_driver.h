//! \file nek_driver.h
//! Driver to initialze and run Nek5000 in stages
#ifndef ENRICO_NEK_DRIVER_H
#define ENRICO_NEK_DRIVER_H

#include "enrico/geom.h"
#include "enrico/heat_fluids_driver.h"
#include "mpi.h"
#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

#include <string>
#include <vector>

namespace enrico {

//! Driver to initialze and run Nek5000 in stages.
class NekDriver : public HeatFluidsDriver {
public:
  //! Initializes Nek5000 with the given MPI communicator.
  //!
  //! A wrapper for the nek_init() routine in Nek5000.  Also retreives common-block
  //! variables from Nek5000 and uses them to set the corresponding attributes in
  //! NekDriver.
  //!
  //! \param comm  The MPI communicator used to initialze Nek5000
  explicit NekDriver(MPI_Comm comm, double pressure_bc, pugi::xml_node xml_root);

  //! Finalizes Nek5000.
  //!
  //! A wrapper for the nek_end() routine in Nek5000.
  ~NekDriver();

  //! Initializes a trivial runtime datafile for Nek5000.
  //!
  //! Nek5000 must read a file with the casename and working directory. It reads this file
  //! instead of reading command-line arguments or otherwise inferring the working
  //! directory.
  void init_session_name();

  //! Runs all timesteps for a heat/fluid solve in Nek5000.
  //!
  //! A wraper for the nek_solve() routine in libnek5000.  This includes the necessary
  //! initialization and finalization for each step.
  void solve_step() final;

  //! Whether the calling rank has access to the full thermal-hydraulic solution field.
  //! Only Nek's master rank has access to the global data; data on other ranks is empty
  bool has_coupling_data() const final { return comm_.rank == 0; }

  xt::xtensor<double, 1> temperature() const final;

  xt::xtensor<double, 1> density() const final;

  xt::xtensor<int, 1> fluid_mask() const final;

  //! Get the coordinate of a local element's centroid.
  //!
  //! The coordinate is dimensionless.  Its units depend on the unit system used that was
  //! used to setup the Nek problem. The user must handle any necessary conversions.
  //!
  //! \param local_elem The local index of the desired element
  //! \return The dimensionless coordinate of the element's centroid
  Position centroid_at(int local_elem) const;

  //! Get the volume of a local element
  //!
  //! The volume is dimensionless.  Its units depend on the unit system used that was used
  //! to setup the Nek problem. The user must handle any necessary conversions.
  //!
  //! \param local_elem The local index of the desired element
  //! \return The dimensionless Volume of the element
  double volume_at(int local_elem) const;

  //! Get the volume-averaged temperature of a local element
  //!
  //! The returned temperature is dimensionless.  Its units depend on the unit system that
  //! was used to setup the Nek5000 problem. The user must handle any necessary
  //! conversions.
  //!
  //! \param local_elem A local element ID
  //! \return The volume-averaged temperature of the element
  double temperature_at(int local_elem) const;

  //! Return true if a local element is in the fluid region
  //! \param local_elem  A local element ID
  //! \return 1 if the local element is in fluid; 0 otherwise
  int in_fluid_at(int local_elem) const;

  //! Set the heat source for a given local element
  //!
  //! The units of heat must match on the unit system that was used to setup the Nek5000
  //! problem (presumably W/cm^3). The caller must handle any necessary conversions.
  //!
  //! \param local_elem A local element ID
  //! \param heat A heat source term
  //! \return Error code
  int set_heat_source_at(int local_elem, double heat);

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
