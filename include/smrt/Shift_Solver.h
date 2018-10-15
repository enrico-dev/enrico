#ifndef Shift_Solver_h
#define Shift_Solver_h

#include <memory>
#include <vector>

#include "Geometria/rtk/RTK_Geometry.hh"
#include "Shift/mc_physics/SCE_Physics.hh"

#include "Omnibus/driver/Multiphysics_Driver.hh"
#include "Assembly_Model.h"
#include "Neutronics_Solver.h"

namespace stream
{

//===========================================================================//
/*!
 * \class Shift_Solver
 * \brief Neutronics solver running Shift problem
 */
/*!
 * \example shift/test/tstShift_Solver.cc
 *
 * Test of Shift_Solver.
 */
//===========================================================================//

class Shift_Solver : public Neutronics_Solver
{
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

  private:
    // >>> DATA
    SP_Assembly_Model   d_assembly;
    SP_Geometry         d_geometry;
    SP_Omn_Driver       d_driver;

    std::vector<double> d_z_edges;
    double      d_power_norm;
    std::string d_power_tally_name;

  public:

    // Constructor
    Shift_Solver(SP_Assembly_Model          assembly,
                 std::string                shift_input,
                 const std::vector<double>& z_edges);

    void solve(const std::vector<double>& fuel_temperature,
               const std::vector<double>& coolant_density,
                     std::vector<double>& power) override;
private:

    // Add power tally to parameter list
    void add_power_tally(RCP_PL&                    pl,
                         const std::vector<double>& z_edges);
};

//---------------------------------------------------------------------------//
} // end namespace stream

//---------------------------------------------------------------------------//
#endif // Shift_Solver_h

//---------------------------------------------------------------------------//
// end of Shift_Solver.h
//---------------------------------------------------------------------------//
