#ifndef SHIFT_NEK_DRIVER_H
#define SHIFT_NEK_DRIVER_H

#include <memory>
#include <vector>

#include "Nemesis/comm/global.hh"

#include "Assembly_Model.h"
#include "Shift_Solver.h"
#include "enrico/error.h"
#include "enrico/message_passing.h"
#include "enrico/nek_driver.h"
#include "nek5000/core/nek_interface.h"

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
class ShiftNekDriver {
private:
  //
  // Data
  //

  // Shift solver
  std::shared_ptr<Shift_Solver> d_shift_solver;

  // Nek solver
  std::shared_ptr<NekDriver> d_nek_solver;

  // Number of local and global T/H elements
  int d_th_num_local;
  int d_th_num_global;

  // MPI datatype for Position objects
  MPI_Datatype d_position_mpi_type;

  // Field data
  std::vector<double> d_temperatures;
  std::vector<double> d_densities;
  std::vector<double> d_powers;

  // Normalization factor for power (average
  double d_power_norm;

public:
  // Constructor
  ShiftNekDriver(std::shared_ptr<Assembly_Model> assembly,
                 const std::vector<double>& z_edges,
                 const std::string& shift_filename,
                 const std::string& enrico_filename,
                 double power_norm,
                 MPI_Comm neutronics_comm,
                 MPI_Comm th_comm);

  // Destructor
  ~ShiftNekDriver();

  // Solve coupled problem
  void solve();

private:
  // Normalize power distribution to appropriate value
  void normalize_power();

  // Map "standard" types to corresponding MPI types
  template<typename T>
  MPI_Datatype get_mpi_type() const
  {
    return MPI_DATATYPE_NULL;
  }

  // Set up MPI types (i.e., Position)
  void init_mpi_datatypes();

  // Free MPI types
  void free_mpi_datatypes();

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
