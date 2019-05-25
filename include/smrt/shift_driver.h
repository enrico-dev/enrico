#ifndef SHIFT_DRIVER_H
#define SHIFT_DRIVER_H

#include <memory>
#include <vector>

#include "Geometria/rtk/RTK_Geometry.hh"
#include "Omnibus/driver/Multiphysics_Driver.hh"
#include "Shift/mc_physics/SCE_Physics.hh"

#include "Assembly_Model.h"
#include "Neutronics_Solver.h"
#include "enrico/geom.h"

namespace enrico {

//===========================================================================//
/*!
 * \class ShiftDriver
 * \brief Neutronics solver running Shift problem
 */
/*!
 * \example shift/test/tstShift_Solver.cc
 *
 * Test of ShiftDriver.
 */
//===========================================================================//

class ShiftDriver : public Neutronics_Solver {
public:
  //@{
  //! Public type aliases
  using SP_Assembly_Model = std::shared_ptr<Assembly_Model>;
  using Omn_Driver = omnibus::Multiphysics_Driver;
  using SP_Omn_Driver = std::shared_ptr<Omn_Driver>;
  using Geometry = geometria::RTK_Core;
  using SP_Geometry = std::shared_ptr<Geometry>;
  using SP_Composition = std::shared_ptr<robus::Composition>;
  using Vec_Composition = std::vector<SP_Composition>;
  using RCP_PL = Teuchos::RCP<Teuchos::ParameterList>;
  //@}

  void update_temperature(const std::vector<double>& temperatures) override;

private:
  // >>> DATA
  SP_Assembly_Model d_assembly;
  SP_Geometry d_geometry;
  SP_Omn_Driver d_driver;

  std::vector<double> d_z_edges;
  std::string d_power_tally_name;

  // Matids corresponding to T/H mesh elements
  int d_num_materials;
  std::vector<int> d_matids;

  // Map from Shift geometric cells to T/H elements
  int d_num_shift_cells;
  std::vector<std::vector<int>> d_power_map;

  // Volume fraction for each T/H element for normalization
  std::vector<double> d_vfracs;

public:
  // Constructor
  ShiftDriver(SP_Assembly_Model assembly,
               std::string shift_input,
               const std::vector<double>& z_edges);

  // Locate centroids from fluids problem
  void set_centroids_and_volumes(const std::vector<enrico::Position>& centroids,
                                 const std::vector<double>& volumes);

  // Solve for new power distribution given temperatures and densities
  void solve(const std::vector<double>& th_temperature,
             const std::vector<double>& coolant_density,
             std::vector<double>& power) override;

private:
  // Add power tally to parameter list
  void add_power_tally(RCP_PL& pl, const std::vector<double>& z_edges);
};

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
#endif // SHIFT_DRIVER_H

//---------------------------------------------------------------------------//
// end of shift_driver.h
//---------------------------------------------------------------------------//
