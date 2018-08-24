#ifndef STREAM_NEK_DRIVER_H
#define STREAM_NEK_DRIVER_H

#include "base_drivers.h"
#include "mpi.h"
#include "stream_geom.h"

namespace stream {

//! Driver to initialze and run Nek5000 in stages.
class NekDriver : public HeatFluidsDriver {
public:

  //! Initializes Nek5000 with the given MPI communicator.
  //!
  //! A wrapper for the nek_init() routine in Nek5000.  Also retreives common-block variables
  //! from Nek5000 and uses them to set the corresponding attributes in NekDriver.
  //!
  //! \param comm  The MPI communicator used to initialze Nek5000
  explicit NekDriver(MPI_Comm comm);

  //! Finalizes Nek5000.
  //!
  //! A wrapper for the nek_end() routine in Nek5000.
  ~NekDriver();

  //! Does nothing; its functionality is not necessary for NekDriver.
  //!
  //! This must be implemented since the virtual method HeatFluidsDriver::init_step() is declared
  //! in the base class.  However, all the actual work for initializing a solve is done in
  //! NekDriver::solve_step().
  void init_step();

  //! Runs all timesteps for a heat/fluid solve in Nek5000.
  //!
  //! A wraper for the nek_solve() routine in libnek5000.  This includes the necessary
  //! initialization/finalization, so NekDriver::init_step() and NekDriver::solve_step() need not
  //! do anything.
  void solve_step();

  //! Does nothing; its functionality is not necessary for NekDriver.
  //!
  //! This must be implemented since the virtual method HeatFluidsDriver::init_step() is declared
  //! in the base class.  However, all the actual work for finalizing a solve is done in
  //! NekDriver::solve_step().
  void finalize_step();

  //! Get the coordinate of a global element's centroid.
  //!
  //! The coordinate is dimensionless.  Its units depend on the unit system used that was used to
  //! setup the Nek problem. The user must handle any necessary conversions.
  //!
  //! \param global_elem The global index of the desired element
  //! \return The dimensionless coordinate of the element's centroid
  Position get_global_elem_centroid(int global_elem) const;

  int lelg_; //!< upper bound on number of mesh elements
  int lelt_; //!< upper bound on number of mesh elements per rank
  int lx1_; //!< polynomial order of the solution
  int nelgt_; //!< total number of mesh elements
  int nelt_; //!< number of local mesh elements
};

} // namespace stream

#endif //STREAM_NEK_DRIVER_H
