#ifndef SHIFT_NEK_DRIVER_H
#define SHIFT_NEK_DRIVER_H

#include <memory>
#include <vector>

#include "Nemesis/comm/global.hh"

#include "Assembly_Model.h"
#include "enrico/error.h"
#include "enrico/message_passing.h"
#include "enrico/nek_driver.h"
#include "nek5000/core/nek_interface.h"
#include "shift_driver.h"
#include "smrt_coupled_driver.h"

namespace enrico {
//===========================================================================//
/*!
 * \class ShiftNekDriver
 * \brief Class for coupling Shift and Nek.
 *
 * This class will perform dampled Picard iteration to converge the
 * coupled nonlinear system.
 */
//===========================================================================//
class ShiftNekDriver : public SmrtCoupledDriver {
public:
  //! Power in [W]
  double power_;

  //! Maximum number of Picard iterations
  int max_picard_iter_;

private:
  //
  // Data
  //

  // Shift solver
  std::shared_ptr<ShiftDriver> d_shift_solver;

  // Nek solver
  std::shared_ptr<NekDriver> d_nek_solver;

  // Number of local and global T/H elements
  int d_th_num_local;
  int d_th_num_global;

  // Field data
  std::vector<double> d_temperatures;
  std::vector<double> d_densities;
  std::vector<double> d_powers;

  // Power indexed by shift cell ID
  std::vector<double> d_power_shift;

public:
  // Constructor
  ShiftNekDriver(std::shared_ptr<Assembly_Model> assembly,
                 const std::vector<double>& z_edges,
                 const std::string& shift_filename,
                 MPI_Comm neutronics_comm,
                 MPI_Comm th_comm);

  // Destructor
  ~ShiftNekDriver();

  // Solve coupled problem
  void solve();

private:
  //
  // T/H communication operations -- ultimately we want to avoid storing
  // the entire T/H solution on a single rank
  //

  // Gather local distributed field into global replicated field
  template<typename T>
  std::vector<T> local_to_global(const std::vector<T>& local_field) const;

  // Scatter global replicated field into local distributed field
  template<typename T>
  std::vector<T> global_to_local(const std::vector<T>& global_field) const;
};

} // end namespace enrico

#endif // SHIFT_NEK_DRIVER_H
