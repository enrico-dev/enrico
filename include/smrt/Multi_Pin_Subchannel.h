#ifndef Multi_Pin_Subchannel_h
#define Multi_Pin_Subchannel_h

#include <memory>
#include <vector>

// Trilinos includes
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

// enrico includes
#include "Assembly_Model.h"
#include "Single_Pin_Subchannel.h"

namespace enrico {

//===========================================================================//
/*!
 * \class Multi_Pin_Subchannel
 * \brief Solver channel-centered subchannel equations over multiple pins.
 */
//===========================================================================//

class Multi_Pin_Subchannel {
public:
  //@{
  //! Typedefs
  using SP_Assembly = std::shared_ptr<Assembly_Model>;
  using RCP_PL = Teuchos::RCP<Teuchos::ParameterList>;
  using SP_Pin_Subchannel = std::shared_ptr<Single_Pin_Subchannel>;
  //@}

private:
  // >>> DATA
  SP_Assembly d_assembly;

  // Note that for channel-centered equations, Nx and Ny are one greater
  //  than number of pins in each direction.
  int d_Nx, d_Ny, d_Nz;

  // Cross sectional area of each channel
  std::vector<double> d_areas;

  // Mass flow rate (kg/s) in each channel
  std::vector<double> d_mdots;

  SP_Pin_Subchannel d_pin_subchannel;

public:
  // Constructor
  Multi_Pin_Subchannel(SP_Assembly assembly,
                       RCP_PL params,
                       const std::vector<double>& dz);

  // Solve
  void solve(const std::vector<double>& power,
             std::vector<double>& channel_temp,
             std::vector<double>& channel_density);

private:
  int channel_index(int ix, int iy) const
  {
    Require(ix < d_Nx);
    Require(iy < d_Ny);
    return ix + d_Nx * iy;
  }
};

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
#endif // Multi_Pin_Subchannel_h

//---------------------------------------------------------------------------//
// end of Multi_Pin_Subchannel.h
//---------------------------------------------------------------------------//
