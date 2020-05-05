//! \file coupled_driver.h
//! Base class for driver that controls a coupled physics solve involving neutronics
//! and thermal-hydraulics physics.
#ifndef ENRICO_COUPLED_DRIVER_H
#define ENRICO_COUPLED_DRIVER_H

#include "enrico/driver.h"
#include "enrico/heat_fluids_driver.h"
#include "enrico/neutronics_driver.h"

#include <pugixml.hpp>
#include <xtensor/xtensor.hpp>

#include <memory> // for unique_ptr
#include <unordered_map>
#include <vector>

namespace enrico {

//! Base class for driver that controls a coupled physics solve involving neutronics
//! and thermal-hydraulics physics.
class CoupledDriver {
public:
  // Types, aliases
  enum class Norm { L1, L2, LINF }; //! Types of norms

  //! Enumeration of available temperature initial condition specifications.
  //! 'neutronics' sets temperature condition from the neutronics input files,
  //! while 'heat' sets temperature based on a thermal-fluids input (or restart) file.
  enum class Initial { neutronics, heat };

  //! Initializes coupled neutron transport and thermal-hydraulics solver with
  //! the given MPI communicator
  //!
  //! \param comm The MPI communicator used for the coupled driver
  //! \param node XML node containing settings
  explicit CoupledDriver(MPI_Comm comm, pugi::xml_node node);

  ~CoupledDriver() {}

  //! Execute the coupled driver
  virtual void execute();

  //! Update the heat source for the thermal-hydraulics solver
  //!
  //! \param relax Apply relaxation to heat source before updating heat solver
  void update_heat_source(bool relax);

  //! Update the temperature for the neutronics solver
  //!
  //! \param relax Apply relaxation to temperature before updating neutronics solver
  void update_temperature(bool relax);

  //! Update the density for the neutronics solver
  //!
  //! \param relax Apply relaxation to density before updating neutronics solver
  void update_density(bool relax);

  //! Check convergence of the coupled solve for the current Picard iteration.
  bool is_converged();

  //! Compute the norm of the temperature between two successive Picard iterations
  //! \param norm enumeration of norm to compute
  //! \return norm of the temperature between two iterations
  double temperature_norm(Norm n);

  //! Get reference to neutronics driver
  //! \return reference to driver
  NeutronicsDriver& get_neutronics_driver() const { return *neutronics_driver_; }

  //! Get reference to thermal-fluids driver
  //! \return reference to driver
  HeatFluidsDriver& get_heat_driver() const { return *heat_fluids_driver_; }

  //! Get timestep iteration index
  //! \return timestep iteration index
  int get_timestep_index() const { return i_timestep_; }

  //! Get Picard iteration index within current timestep
  //! \return Picard iteration index within current timestep
  int get_picard_index() const { return i_picard_; }

  //! Whether solve is for first Picard iteration of first timestep
  bool is_first_iteration() const
  {
    return get_timestep_index() == 0 and get_picard_index() == 0;
  }

  Comm comm_; //!< The MPI communicator used to run the driver

  double power_; //!< Power in [W]

  int max_timesteps_; //!< Maximum number of time steps

  int max_picard_iter_; //!< Maximum number of Picard iterations

  //! Picard iteration convergence tolerance, defaults to 1e-3 if not set
  double epsilon_{1e-3};

  //! Constant relaxation factor for the heat source,
  //! defaults to 1.0 (standard Picard) if not set
  double alpha_{1.0};

  //! Constant relaxation factor for the temperature, defaults to the
  //! relaxation aplied to the heat source if not set
  double alpha_T_{alpha_};

  //! Constant relaxation factor for the density, defaults to the
  //! relaxation applied to the heat source if not set
  double alpha_rho_{alpha_};

  //! Where to obtain the temperature initial condition from. Defaults to the
  //! temperatures in the neutronics input file.
  Initial temperature_ic_{Initial::neutronics};

  //! Where to obtain the density initial condition from. Defaults to the densities
  //! in the neutronics input file.
  Initial density_ic_{Initial::neutronics};

private:
  //! Create bidirectional mappings from neutronics cell instances to/from TH elements
  void init_mappings();

  //! Initialize the Monte Carlo tallies for all cells
  void init_tallies();

  //! Initialize global volume buffers for neutronics ranks
  void init_volumes();

  //! Initialize global fluid masks on all TH ranks.
  void init_elem_fluid_mask();

  //! Initialize fluid masks for neutronics cells on all neutronic ranks.
  void init_cell_fluid_mask();

  //! Initialize current and previous Picard temperature fields
  void init_temperatures();

  //! Initialize current and previous Picard density fields
  void init_densities();

  //! Initialize current and previous Picard heat source fields. Note that
  //! because the neutronics solver is assumed to run first, that no initial
  //! condition is required for the heat source. So, unlike init_temperatures(),
  //! this method does not set any initial values.
  void init_heat_source();

  //! Print report of communicator layout
  void comm_report();

  //! Special alpha value indicating use of Robbins-Monro relaxation
  constexpr static double ROBBINS_MONRO = -1.0;

  int i_timestep_; //!< Index pertaining to current timestep

  int i_picard_; //!< Index pertaining to current Picard iteration

  //! The rank in comm_ that corresponds to the root of the neutronics comm
  int neutronics_root_ = MPI_PROC_NULL;

  //! The rank in comm_ that corresponds to the root of the heat comm
  int heat_root_ = MPI_PROC_NULL;

  //! Current Picard iteration temperature; this temperature is the temperature
  //! computed by the thermal-hydraulic solver, and data mappings may result in
  //! a different temperature actually used in the neutronics solver. For example,
  //! the entries in this xtensor may be averaged over neutronics cells to give
  //! the temperature used by the neutronics solver.
  xt::xtensor<double, 1> temperatures_;

  xt::xtensor<double, 1> temperatures_prev_; //!< Previous Picard iteration temperature

  //! Current Picard iteration density; this density is the density
  //! computed by the thermal-hydraulic solver, and data mappings may result in
  //! a different density actually used in the neutronics solver. For example,
  //! the entries in this xtensor may be averaged over neutronics cells to give
  //! the density used by the neutronics solver.
  xt::xtensor<double, 1> densities_;

  xt::xtensor<double, 1> densities_prev_; //!< Previous Picard iteration density

  //! Current Picard iteration heat source; this heat source is the heat source
  //! computed by the neutronics solver, and data mappings may result in a different
  //! heat source actually used in the heat solver. For example, the entries in this
  //! xtensor may be averaged over thermal-hydraulics cells to give the heat source
  //! used by the thermal-hydraulics solver.
  xt::xtensor<double, 1> heat_source_;

  xt::xtensor<double, 1> heat_source_prev_; //!< Previous Picard iteration heat source

  std::unique_ptr<NeutronicsDriver> neutronics_driver_;  //!< The neutronics driver
  std::unique_ptr<HeatFluidsDriver> heat_fluids_driver_; //!< The heat-fluids driver

  //! States whether a global element is in the fluid region
  //! These are **not** ordered by TH global element indices.  Rather, these are
  //! ordered according to an MPI_Gatherv operation on TH local elements.
  std::vector<int> elem_fluid_mask_;

  //! States whether a neutronic cell is in the fluid region
  xt::xtensor<int, 1> cell_fluid_mask_;

  //! Volumes of global elements in TH solver
  //! These are **not** ordered by TH global element indices.  Rather, these are
  //! ordered according to an MPI_Gatherv operation on TH local elements.
  std::vector<double> elem_volumes_;

  //! Map that gives a list of TH element indices for a given neutronics cell
  //! handle. The TH element indices refer to indices defined by the MPI_Gatherv
  //! operation, and do not reflect TH internal global element indexing.
  std::unordered_map<CellHandle, std::vector<int32_t>> cell_to_elems_;

  //! Map that gives the neutronics cell handle for a given TH element index.
  //! The TH element indices refer to indices defined by the MPI_Gatherv
  //! operation, and do not reflect TH internal global element indexing.
  std::vector<CellHandle> elem_to_cell_;

  //! Number of cell instances in neutronics model
  int32_t n_cells_;

  //! Number of global elements in heat/fluids model
  int32_t n_global_elem_;

  // Norm to use for convergence checks
  Norm norm_{Norm::LINF};
};

} // namespace enrico

#endif // ENRICO_COUPLED_DRIVER_H
