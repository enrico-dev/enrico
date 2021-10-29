//! \file coupled_driver.h
//! Base class for driver that controls a coupled physics solve involving neutronics
//! and thermal-hydraulics physics.
#ifndef ENRICO_COUPLED_DRIVER_H
#define ENRICO_COUPLED_DRIVER_H

#include "enrico/driver.h"
#include "enrico/heat_fluids_driver.h"
#include "enrico/neutronics_driver.h"
#include "enrico/timer.h"

#include <pugixml.hpp>
#include <xtensor/xtensor.hpp>

#include <map>
#include <memory> // for unique_ptr
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
  CoupledDriver(MPI_Comm comm, pugi::xml_node node);

  ~CoupledDriver() = default;

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

  //! Report cumulative times for CoupledDriver member functions
  void timer_report();

  Timer timer_init_comms;         //!< For initialzing subcommunicators, etc.
  Timer timer_init_mapping;       //!< For the init_mapping() member function
  Timer timer_init_tallies;       //!< For the init_tallies() member function
  Timer timer_init_volume;        //!< For the init_volume() member function
  Timer timer_init_fluid_mask;    //!< For the init_fluid_mask() member function
  Timer timer_init_temperature;   //!< For the init_temperature() member function
  Timer timer_init_density;       //!< For the init_density() member function
  Timer timer_init_heat_source;   //!< For the init_heat_source() member function
  Timer timer_update_density;     //!< For the update_density() member function
  Timer timer_update_heat_source; //!< For the update_heat_source() member function
  Timer timer_update_temperature; //!< For the update_temperature() member function

private:
  //! Parse coupled driver's runtime parameters from enrico.xml
  void parse_xml_params(const pugi::xml_node& node);

  //! Create subcommunicators for single-physics drivers
  void init_comms(const pugi::xml_node& node);

  //! Create mappings between neutronics cell instances and heat/fluids elements
  void init_mapping();

  //! Initialize the Monte Carlo tallies for all cells
  void init_tallies();

  //! Calculate and store local cell volumes in each heat/fluids rank
  void init_volume();

  //! Report how closely the neutron driver's volumes and the calculated local cell
  //! volumes (from init_volume()) match.  Raises no errors or warnings.
  void check_volumes();

  //! Initialize fluid mask for local cells on each heat/fluids rank
  void init_fluid_mask();

  //! Initialize current and previous Picard temperature fields
  void init_temperature();

  //! Initialize current and previous Picard density fields
  void init_density();

  //! Initialize current and previous Picard heat source fields.
  //!
  //! Because the neutronics solver is assumed to run first, no initial
  //! condition is required for the heat source. Thus, unlike init_temperature(),
  //! this member function does not set any initial values.
  void init_heat_source();

  //! Print report of communicator layout if high verbosity is set
  void comm_report();

  //! Special alpha value indicating use of Robbins-Monro relaxation
  constexpr static double ROBBINS_MONRO = -1.0;

  int i_timestep_; //!< Index pertaining to current timestep

  int i_picard_; //!< Index pertaining to current Picard iteration

  //! The rank in comm_ that corresponds to the root of the neutronics comm
  int neutronics_root_ = MPI_PROC_NULL;

  //! The rank in comm_ that corresponds to the root of the heat comm
  int heat_root_ = MPI_PROC_NULL;

  //! List of ranks in this->comm_ that are in the heat/fluids subcomm
  std::vector<int> heat_ranks_;

  //! List of ranks in this->comm_ that are in the neutronics subcomm
  std::vector<int> neutronics_ranks_;

  //! Local cell temperature at current Picard iteration. Set only on heat/fluids ranks.
  xt::xtensor<double, 1> cell_temperature_;

  //! Local cell temperature at previous Picard iteration. Set only on heat/fluids ranks.
  xt::xtensor<double, 1> cell_temperature_prev_;

  //! Local cell density at current Picard iteration. Set only on heat/fluids ranks.
  xt::xtensor<double, 1> cell_density_;

  //! Local cell density at previous Picard iteration. Set only on heat/fluids ranks.
  xt::xtensor<double, 1> cell_density_prev_;

  //! Local cell heat source at current Picard iteration. Set only on heat/fluids ranks.
  xt::xtensor<double, 1> cell_heat_source_;

  //! Local cell heat source at previous Picard iteration. Set only on heat/fluids ranks.
  xt::xtensor<double, 1> cell_heat_source_prev_;

  std::unique_ptr<NeutronicsDriver> neutronics_driver_;  //!< The neutronics driver
  std::unique_ptr<HeatFluidsDriver> heat_fluids_driver_; //!< The heat-fluids driver

  //! 1 if local cell is in fluid, 0 if in solid. Set only on heat/fluids ranks.
  std::vector<int> cell_fluid_mask_;

  //! Maps heat/fluids element index to global cell handle.  Set only on heat/fluids ranks.
  std::vector<CellHandle> elem_to_glob_cell_;

  //! Maps local cell index to global cell handle.  Set only on heat/fluid ranks.
  std::vector<CellHandle> cell_to_glob_cell_;

  //! Maps global cell handle to local element index.  Set only on heat/fluids ranks.
  std::map<CellHandle, std::vector<int32_t>> glob_cell_to_elem_;

  //! Local cell volumes.  Set only on heat/fluids ranks.
  std::vector<double> cell_volume_;

  //! Local element volumes.  Set only on heat/fluids ranks.
  std::vector<double> elem_volume_;

  // Norm to use for convergence checks
  Norm norm_{Norm::LINF};

  // Print verbose output
  bool verbose_ = false;
};

} // namespace enrico

#endif // ENRICO_COUPLED_DRIVER_H
