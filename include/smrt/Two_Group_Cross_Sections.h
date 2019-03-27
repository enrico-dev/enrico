#ifndef Two_Group_Cross_Sections_h
#define Two_Group_Cross_Sections_h

#include "Assembly_Model.h"

namespace enrico
{

//===========================================================================//
/*!
 * \class Two_Group_Cross_Sections
 * \brief Class for interpolated two group cross section data.
 */
//===========================================================================//

class Two_Group_Cross_Sections
{
  public:

    struct XS_Data
    {
        double diffusion[2];
        double absorption[2];
        double scatter;
        double nu_fission[2];
    };

  private:
    // >>> DATA

  public:

    // Constructor
    Two_Group_Cross_Sections(){}

    // Get XS data at given temperature and density
    XS_Data get_data(Assembly_Model::PIN_TYPE type, double T, double rho);
};

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
#endif // Two_Group_Cross_Sections_h

//---------------------------------------------------------------------------//
// end of Two_Group_Cross_Sections.h
//---------------------------------------------------------------------------//
