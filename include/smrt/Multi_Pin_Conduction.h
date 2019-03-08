#ifndef Multi_Pin_Conduction_h
#define Multi_Pin_Conduction_h

#include <memory>

// Trilinos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// stream includes
#include "Assembly_Model.h"
#include "Single_Pin_Conduction.h"

namespace stream
{

//===========================================================================//
/*!
 * \class Multi_Pin_Conduction
 * \brief Solve heat equation over multiple fuel pins
 */
/*!
 * \example Surrogate/test/tstMulti_Pin_Conduction.cc
 *
 * Test of Multi_Pin_Conduction.
 */
//===========================================================================//

class Multi_Pin_Conduction
{
  public:
    //@{
    //! Typedefs
    using SP_Assembly = std::shared_ptr<Assembly_Model>;
    using SP_Pin_Conduction = std::shared_ptr<Single_Pin_Conduction>;
    using RCP_PL = Teuchos::RCP<Teuchos::ParameterList>;
    //@}

  private:
    // >>> DATA
    SP_Assembly d_assembly;
    int d_Nz;

    // Single pin solver
    SP_Pin_Conduction d_pin_conduction;

  public:

    // Constructor
    Multi_Pin_Conduction(SP_Assembly                assembly,
                         RCP_PL                     parameters,
                         const std::vector<double>& dz);

    // Solve for all pins
    void solve(const std::vector<double>& power,
               const std::vector<double>& channel_temp,
                     std::vector<double>& fuel_temp);
};

//---------------------------------------------------------------------------//
} // end namespace stream

//---------------------------------------------------------------------------//
#endif // Multi_Pin_Conduction_h

//---------------------------------------------------------------------------//
// end of Multi_Pin_Conduction.h
//---------------------------------------------------------------------------//
