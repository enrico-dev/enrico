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
  virtual void update_heat_source() final;

  //! Set the heat source in the thermal-hydraulics solver
  virtual void set_heat_source() {}

  //! Update the temperature for the neutronics solver
  virtual void update_temperature() {}

  //! Update the density for the neutronics solver
  virtual void update_density() {}

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

  //! Constant relaxation factor, defaults to 1.0 (standard Picard) if not set
  double alpha_{1.0};

  //! Enumeration of available temperature initial condition specifications.
  //! 'neutronics' sets temperature condition from the neutronics input files,
  //! while 'heat' sets temperature based on a thermal-fluids input (or restart) file.
  enum class Initial {neutronics, heat};

  //! Where to obtain the temperature initial condition from. Defaults to the
  //! temperatures in the neutronics input file.
  Initial temperature_ic_{Initial::neutronics};

protected:
  //! Initialize current and previous Picard temperature fields
  virtual void init_temperatures() {}

  //! Initialize current and previous Picard heat source fields. Note that
  //! because the neutronics solver is assumed to run first, that no initial
  //! condition is required for the heat source. So, unlike init_temperatures(),
  //! this method does not set any initial values.
  virtual void init_heat_source() {}

  xt::xtensor<double, 1> temperatures_; //!< Current Picard iteration temperature

  xt::xtensor<double, 1> temperatures_prev_; //!< Previous Picard iteration temperature

  xt::xtensor<double, 1> heat_source_; //!< Current Picard iteration heat source

  xt::xtensor<double, 1> heat_source_prev_; //!< Previous Picard iteration heat source

private:
  int i_timestep_; //!< Index pertaining to current timestep

  int i_picard_; //!< Index pertaining to current Picard iteration
};

} // namespace enrico

#endif // ENRICO_COUPLED_DRIVER_H
