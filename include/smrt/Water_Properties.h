#ifndef Water_Properties_h
#define Water_Properties_h

namespace enrico
{

//===========================================================================//
/*!
 * \class Water_Properties
 * \brief Get material properties for water
 */
//===========================================================================//

class Water_Properties
{
  public:

    static double Density(double h, double p);
    static double Enthalpy(double T, double p);
    static double Temperature(double h, double p);
};

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
#endif // Water_Properties_h

//---------------------------------------------------------------------------//
// end of Water_Properties.h
//---------------------------------------------------------------------------//
