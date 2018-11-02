//! \file openmc_nek_driver.h
//! Driver for coupled Nek5000/OpenMC simulations
#ifndef STREAM_OPENMC_NEK_DRIVER_H
#define STREAM_OPENMC_NEK_DRIVER_H

#include "base_drivers.h"
#include "mpi.h"
#include "nek_driver.h"
#include "openmc_driver.h"

#include <memory> // for unique_ptr
#include <unordered_set>

namespace stream {

//! Driver for coupled Nek5000/OpenMC simulations
//!
//! \todo This is not actually derived from CoupledDriver.  Currently, it is unclear how or if the
//! base class will be implemented.  The issue will be revisited
class OpenmcNekDriver {
public:
  //! Given existing MPI comms, initialize drivers and geometry mappings
  //!
  //! Currently, openmc_comm and nek_comm must be subsets of coupled_comm.  The function
  //! stream::get_node_comms() can be used to split a coupled_comm into suitable subcomms.
  //!
  //! \param power Power in [W]
  //! \param coupled_comm An existing communicator for the coupled driver
  OpenmcNekDriver(MPI_Comm coupled_comm, pugi::xml_node xml_root);

  //! Frees any data structures that need manual freeing.
  ~OpenmcNekDriver();

  //! Transfers heat source terms from Nek5000 to OpenMC
  void update_heat_source();

  //! Tranfsers temperatures from OpenMC to Nek5000
  void update_temperature();

  //! Run one timstep
  void solve_in_time();

  Comm comm_; //!< The communicator used to run this driver
  Comm intranode_comm_;  //!< The communicator reprsenting intranode ranks
  std::unique_ptr<OpenmcDriver> openmc_driver_;  //!< The OpenMC driver
  std::unique_ptr<NekDriver> nek_driver_;  //!< The Nek5000 driver
  double power_; //!< Power in [W]
  int max_timesteps_; //! Maximum of timesteps
  int max_picard_iter_; //! Maximum number of Picard iterations per timestep
private:

  //! Initialize MPI datatypes (currently, only position_mpi_datatype)
  void init_mpi_datatypes();

  //! Create bidirectional mappings from OpenMC materials to/from Nek5000 elements
  void init_mappings();

  //! Initialize the tallies for all OpenMC materials
  void init_tallies();

  //! Initialize global temperature buffers for OpenMC ranks
  void init_temperatures();

  //! Initialize global volume buffers for OpenMC ranks
  void init_volumes();

  //! Get the heat index for a given OpenMC material
  //! \param mat_index An OpenMC material index
  //! \return The heat index
  int get_heat_index(int32_t mat_index) const
  {
    return heat_index_.at(mat_index - 1);
  }

  //! Frees the MPI datatypes (currently, only position_mpi_datatype)
  void free_mpi_datatypes();

  //! MPI datatype for sending/receiving Position objects.
  MPI_Datatype position_mpi_datatype;

  //! Gives a Position of a global element's centroid
  //! These are **not** ordered by Nek's global element indices.  Rather, these are ordered
  //! according to an MPI_Gatherv operation on Nek5000's local elements.
  std::vector<Position> global_elem_centroids_;

  //! The dimensionless temperatures of Nek's global elements
  //! These are **not** ordered by Nek's global element indices.  Rather, these are ordered
  //! according to an MPI_Gatherv operation on Nek5000's local elements.
  std::vector<double> global_elem_temperatures_;

  //! The dimensionless volumes of Nek's global elements
  //! These are **not** ordered by Nek's global element indices.  Rather, these are ordered
  //! according to an MPI_Gatherv operation on Nek5000's local elements.
  std::vector<double> global_elem_volumes_;

  //! Map that gives a list of Nek element global indices for a given OpenMC material index
  std::unordered_map<int32_t, std::vector<int>> mat_to_elems_;

  //! Map that gives the OpenMC material index for a given Nek global element index
  std::unordered_map<int, int32_t> elem_to_mat_;

  //! Mapping of material indices (minus 1) to positions in array of heat sources that is used
  //! during update_heat_source
  std::vector<int> heat_index_;

  //! Number of materials in OpenMC model
  int32_t n_materials_;

  //! Number of Nek local elements on this MPI rank.
  //! If nek_driver_ is active, this equals nek_driver.nelt_.  If not, it equals 0.
  int n_local_elem_;

  //! Number of Nek global elements across all ranks.
  //! Always equals nek_driver_.nelgt_.
  int n_global_elem_;

};

} // namespace stream

#endif //STREAM_OPENMC_NEK_DRIVER_H
