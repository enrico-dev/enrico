//! \file coupled_driver.h
//! Base class for driver that controls a coupled physics solve involving neutronics
//! and thermal-hydraulics physics.
#ifndef ENRICO_COUPLED_DRIVER_H
#define ENRICO_COUPLED_DRIVER_H

#include "enrico/driver.h"
#include "enrico/neutronics_driver.h"
#include "heat_fluids_driver.h"
#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

#include <memory>

namespace enrico {

//! Base class for driver that controls a coupled physics solve involving neutronics
//! and thermal-hydraulics physics.
class CoupledDriver {
public:
  //! Initializes coupled neutron transport and thermal-hydraulics solver with
  //! the given MPI communicator
  //!
  //! \param comm The MPI communicator used for the coupled driver
  //! \param node XML node containing settings
  explicit CoupledDriver(MPI_Comm comm, pugi::xml_node node);

  ~CoupledDriver() {}

  //! Execute the coupled driver
  virtual void execute();

  //! Whether the calling rank has access to the coupled solution
  //! fields for the heat source, temperature, density, and other protected member
  //! variables of this class.
  virtual bool has_global_coupling_data() const = 0;

  //! Update the heat source for the thermal-hydraulics solver
  void update_heat_source();

  //! Set the heat source in the thermal-hydraulics solver
  virtual void set_heat_source() {}

  //! Update the temperature for the neutronics solver
  void update_temperature();

  //! Set the temperature in the neutronics solver
  virtual void set_temperature() {};

  //! Update the density for the neutronics solver
  void update_density();

  //! Update the density for the neutronics solver
  virtual void set_density() {}

  //! Check convergence of the coupled solve for the current Picard iteration.
  virtual bool is_converged();

  enum class Norm { L1, L2, LINF }; //! Types of norms

  //! Compute the norm of the temperature between two successive Picard iterations
  //! \param n enumeration of norm to compute
  //! \return norm of the temperature between two iterations
  //! \return whether norm is less than convergence tolerance
  void compute_temperature_norm(const Norm& n, double& norm, bool& converged);

  //! Get reference to neutronics driver
  //! \return reference to driver
  virtual NeutronicsDriver& get_neutronics_driver() const = 0;

  //! Get reference to thermal-fluids driver
  //! \return reference to driver
  virtual HeatFluidsDriver& get_heat_driver() const = 0;

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

  //! Enumeration of available temperature initial condition specifications.
  //! 'neutronics' sets temperature condition from the neutronics input files,
  //! while 'heat' sets temperature based on a thermal-fluids input (or restart) file.
  enum class Initial {neutronics, heat};

  //! Where to obtain the temperature initial condition from. Defaults to the
  //! temperatures in the neutronics input file.
  Initial temperature_ic_{Initial::neutronics};

  //! Where to obtain the density initial condition from. Defaults to the densities
  //! in the neutronics input file.
  Initial density_ic_{Initial::neutronics};

protected:
  //! Initialize current and previous Picard temperature fields
  virtual void init_temperatures() {}

  //! Initialize current and previous Picard density fields
  virtual void init_densities() {}

  //! Initialize current and previous Picard heat source fields. Note that
  //! because the neutronics solver is assumed to run first, that no initial
  //! condition is required for the heat source. So, unlike init_temperatures(),
  //! this method does not set any initial values.
  virtual void init_heat_source() {}

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

private:
  int i_timestep_; //!< Index pertaining to current timestep

  int i_picard_; //!< Index pertaining to current Picard iteration
};

} // namespace enrico

#endif // ENRICO_COUPLED_DRIVER_H
