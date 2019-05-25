//! \file openmc_nek_driver.h
//! Driver for coupled Nek5000/OpenMC simulations
#ifndef ENRICO_OPENMC_NEK_DRIVER_H
#define ENRICO_OPENMC_NEK_DRIVER_H

#include "enrico/coupled_driver.h"
#include "enrico/nek_driver.h"
#include "enrico/openmc_driver.h"
#include "enrico/message_passing.h"
#include "mpi.h"

#include <unordered_set>

namespace enrico {

//! Driver for coupled Nek5000/OpenMC simulations
class OpenmcNekDriver : public CoupledDriver {

public:
  //! Initializes coupled OpenMC-Nek5000 driver with the given MPI communicator and
  //! sets up geometry mappings.
  //!
  //! Currently, openmc_comm and nek_comm must be subsets of comm.  The function
  //! enrico::get_node_comms() can be used to split a coupled_comm into suitable subcomms.
  //!
  //! \param comm The MPI communicator used for the coupled driver
  //! \param node XML node containing settings
  explicit OpenmcNekDriver(MPI_Comm comm, pugi::xml_node node);

  //! Frees any data structures that need manual freeing.
  ~OpenmcNekDriver();

  //! Whether the calling rank has access to global coupling fields. Because the OpenMC
  //! and Nek communicators are assumed to overlap (though they are not the same), and
  //! Nek broadcasts its solution onto the OpenMC ranks, we need to check that both
  //! communicators are active.
  //!
  //! TODO: This won't work if the OpenMC and Nek communicators are disjoint
  bool has_global_coupling_data() const override;

  void set_heat_source() override;

  void update_temperature() override;

  void update_density() override;

  NeutronicsDriver& get_neutronics_driver() const override;

  HeatFluidsDriver & get_heat_driver() const override;

  Comm intranode_comm_; //!< The communicator representing intranode ranks
  std::unique_ptr<OpenmcDriver> openmc_driver_; //!< The OpenMC driver
  std::unique_ptr<NekDriver> nek_driver_;       //!< The Nek5000 driver
  int openmc_procs_per_node_; //!< Number of MPI ranks per (shared-memory) node in OpenMC comm

protected:
  //! Initialize global temperature buffers on all OpenMC ranks.
  //!
  //! These arrays store the dimensionless temperatures of Nek's global elements. These are **not**
  //! ordered by Nek's global element indices. Rather, these are ordered according
  //! to an MPI_Gatherv operation on Nek5000's local elements.
  void init_temperatures() override;

  //! Initialize global source buffers on all OpenMC ranks.
  //!
  //! These arrays store the dimensionless source of Nek's global elements. These are **not**
  //! ordered by Nek's global element indices. Rather, these are ordered according
  //! to an MPI_Gatherv operation on Nek5000's local elements.
  void init_heat_source() override;

  //! Initialize global fluid masks on all OpenMC ranks.
  //!
  //! These arrays store the dimensionless source of Nek's global elements. These are **not**
  //! ordered by Nek's global element indices. Rather, these are ordered according
  //! to an MPI_Gatherv operation on Nek5000's local elements.
  void init_elem_fluid_mask();

  //! Initialize fluid masks for OpenMC cells on all OpenMC ranks.
  void init_cell_fluid_mask();

private:
  //! Initialize MPI datatypes (currently, only position_mpi_datatype)
  void init_mpi_datatypes();

  //! Create bidirectional mappings from OpenMC materials to/from Nek5000 elements
  void init_mappings();

  //! Initialize the tallies for all OpenMC materials
  void init_tallies();

  //! Initialize global volume buffers for OpenMC ranks
  void init_volumes();

  //! Allocate space for the global volume buffers in OpenMC ranks
  //! Currently, the density values are uninitialized.
  void init_elem_densities();

  //! Frees the MPI datatypes (currently, only position_mpi_datatype)
  void free_mpi_datatypes();

  //! Updates elem_densities_ from Nek5000's density data
  void update_elem_densities();

  //! Updates cell_densities_ from elem_densities_.
  void update_cell_densities();

  //! MPI datatype for sending/receiving Position objects.
  MPI_Datatype position_mpi_datatype;

  //! Gives a Position of a global element's centroid
  //! These are **not** ordered by Nek's global element indices.  Rather, these are
  //! ordered according to an MPI_Gatherv operation on Nek5000's local elements.
  std::vector<Position> elem_centroids_;

  //! States whether a global element is in the fluid region
  //! These are **not** ordered by Nek's global element indices.  Rather, these are
  //! ordered according to an MPI_Gatherv operation on Nek5000's local elements.
  xt::xtensor<int, 1> elem_fluid_mask_;

  //! States whether an OpenMC cell in the fluid region
  xt::xtensor<int, 1> cell_fluid_mask_;

  //! The dimensionless volumes of Nek's global elements
  //! These are **not** ordered by Nek's global element indices.  Rather, these are
  //! ordered according to an MPI_Gatherv operation on Nek5000's local elements.
  xt::xtensor<double, 1> elem_volumes_;

  //! The dimensionless densities of Nek's global elements
  //! These are **not** ordered by Nek's global element indices.  Rather, these are
  //! ordered according to an MPI_Gatherv operation on Nek5000's local elements.
  xt::xtensor<double, 1> elem_densities_;

  //! Map that gives a list of Nek element global indices for a given OpenMC material
  //! index. The Nek global element indices refer to indices defined by the MPI_Gatherv
  //! operation, and do not reflect Nek's internal global element indexing.
  std::unordered_map<int32_t, std::vector<int>> mat_to_elems_;

  //! Map that gives the OpenMC material index for a given Nek global element index.
  //! The Nek global element indices refer to indices defined by the MPI_Gatherv
  //! operation, and do not reflect Nek's internal global element indexing.
  std::unordered_map<int, int32_t> elem_to_mat_;

  //! Mapping of material indices (minus 1) to positions in array of heat sources that is
  //! used during update_heat_source
  std::vector<int> heat_index_;

  //! Number of materials in OpenMC model
  int32_t n_materials_;

  //! Number of Nek local elements on this MPI rank.
  //! If nek_driver_ is active, this equals nek_driver.nelt_.  If not, it equals 0.
  size_t n_local_elem_;

  //! Number of Nek global elements across all ranks.
  //! Always equals nek_driver_.nelgt_.
  size_t n_global_elem_;
};

} // namespace enrico

#endif // ENRICO_OPENMC_NEK_DRIVER_H
