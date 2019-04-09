//! \file coupled_driver.h
//! Base class for driver that controls a coupled physics solve involving neutronics
//! and thermal-hydraulics physics.
#ifndef ENRICO_COUPLED_DRIVER_H
#define ENRICO_COUPLED_DRIVER_H

#include "pugixml.hpp"
#include "enrico/driver.h"

#include <memory> // for unique_ptr

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

  ~CoupledDriver() {};

  //! Execute the coupled driver
  virtual void execute();

  //! Update the heat source for the thermal-hydraulics solver
  virtual void update_heat_source() {};

  //! Update the temperature for the neutronics solver
  virtual void update_temperature() {};

  //! Update the density for the neutronics solver
  virtual void update_density() {};

  //! Get reference to neutronics driver
  //! \return reference to driver
  virtual Driver& getNeutronicsDriver() const = 0;

  //! Get reference to thermal-fluids driver
  //! \return reference to driver
  virtual Driver & getHeatDriver() const = 0;

  Comm comm_; //! The MPI communicator used to run the driver

  double power_; //! Power in [W]

  int max_timesteps_; //! Maximum number of time steps

  int max_picard_iter_; //! Maximum number of Picard iterations
};

} // namespace enrico

#endif //ENRICO_COUPLED_DRIVER_H
