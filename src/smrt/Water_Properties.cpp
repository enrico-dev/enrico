#include <cmath>

#include "smrt/Water_Properties.h"

#include "Nemesis/harness/DBC.hh"

namespace enrico
{

namespace
{

double PA_TO_PSI = 0.000145038;
double JG_TO_BTULBM = 0.00043021;
double LBM_TO_KG = 0.453592;
double FT_TO_M = 0.3048;
double PCRIT = 3208.2;
double HCRIT = 906.0;

}

//---------------------------------------------------------------------------//
/* Computes density of liquid water given known specific enthalpy, h,
 * and pressure, p.  Specific enthalpy should be provided in J/kg and
 * pressure in Pa.  Density is returned in kg/m^3.
 * Uses equation III.1-8a and III.1-8b with data from Table III.1-5 of
 * "M. P. Paulsen, et. al, RETRAN-3D,
 * Electric Power Research Institute, Technical Document NP-7450,
 * Volume 1, September 1998"
 * Only the equations for liquid water are implemented and no check is currently
 * performed to ensure that the enthalpy is below the saturation enthalpy.
 */
//---------------------------------------------------------------------------//
double Water_Properties::Density(double h, double p)
{
    double CN0[3][3] =
        {{1.1956901695e-9,   3.7591804374e-11, -2.4473796276e-13},
         {1.6173258743e-13, -2.1383283863e-14,  9.3054844544e-17},
         {7.4927085737e-17,  4.2176141427e-18, -1.1512516405e-20}};

    double CN1[3][5] =
        {{-0.4117961750e1,  -0.3811294543e-3,  0.4308265942e-5,
          -0.9160120130e-8,  0.8017924673e-11},
         {-0.4816067020e-5,  0.7744786733e-7, -0.6988467605e-9,
           0.1916720525e-11,-0.1760288590e-14},
         {-0.1820625039e-8,  0.1440785930e-10,-0.2082170753e-13,
          -0.3603625114e-16, 0.7407124321e-19}};

    // Convert to British
    p *= PA_TO_PSI;
    h *= JG_TO_BTULBM;

    Check(p > 0.01);
    Check(p < PCRIT);
    Check(h > 0.0);
    Check(h < HCRIT);

    double nu = 0.0;
    if (h < 250.0)
    {
        // Equation III.1-8a
        double exp_term = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                nu = nu + CN0[i][j] * std::pow(p,i) * std::pow(250.0 - h,j+2);
            }
        }
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 5; ++j)
            {
                exp_term = exp_term + CN1[i][j] * std::pow(p,i) * std::pow(h,j);
            }
        }
        nu = nu + std::exp(exp_term);
    }
    else
    {
        // Equation III.1-8b
        double exp_term = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 5; ++j)
            {
                exp_term = exp_term + CN1[i][j] * std::pow(p,i) * std::pow(h,j);
            }
        }
        nu = std::exp(exp_term);
    }

    // Convert back to metric
    nu *= std::pow(FT_TO_M,3) / LBM_TO_KG;
    return 1.0 / nu;
}

//---------------------------------------------------------------------------//
/* Computes temperature of liquid water given known specific enthalpy, h,
 * and pressure, p.  Specific enthalpy should be provided in J/kg and
 * pressure in Pa.  Temperature is returned in K.
 * Uses equation III.1-6a with data from Table III.1-4 of
 * "M. P. Paulsen, et. al, RETRAN-3D,
 * Electric Power Research Institute, Technical Document NP-7450,
 * Volume 1, September 1998"
 * Only the equations for liquid water are implemented and no check is currently
 * performed to ensure that the enthalpy is below the saturation enthalpy.
 */
//---------------------------------------------------------------------------//
double Water_Properties::Temperature(double h, double p)
{
    double CT1[2][4] =
       {{0.3276275552e2,   0.9763617000e0,  0.1857226027e-3, -0.4682674330e-6},
        {0.3360880214e-2, -0.5595281760e-4, 0.1618595991e-6, -0.1180204381e-9}};

    // Equations require pressure in psia and specific enthalpy in Btu/lbm
    // and produce a result in degrees F.
    p = p * PA_TO_PSI;
    h = h * JG_TO_BTULBM;

    // Sanity check on inputs
    Check(p > 0.01);
    Check(p < PCRIT);
    Check(h > 0.0);
    Check(h < HCRIT);

    // Evaluate based on formula
    double T = 0.0;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            T = T + CT1[i][j] * std::pow(p,i) * std::pow(h,j);
        }
    }

    // Convert from F to K
    T = (T + 459.67) * 5.0/9.0;
    return T;
}

//---------------------------------------------------------------------------//
/* Computes specific enthalpy of liquid water given known Temperature, T,
 * and pressure, p.  Temperature should be provided in K and pressure in
 * Pa.  Specific enthalpy will be returned in J/kg.
 * Uses equation III.1-6a with data from Table III.1-4 of
 * "M. P. Paulsen, et. al, RETRAN-3D,
 * Electric Power Research Institute, Technical Document NP-7450,
 * Volume 1, September 1998"
 * The above citation only converts from (h,p) to T, and not from (T,p) to h.
 * We use bisection to compute the h from (T,p) using the WaterTemperature
 * function.
 */
//---------------------------------------------------------------------------//
double Water_Properties::Enthalpy(double T, double p)
{
    double hmin = 1.0e-6 / JG_TO_BTULBM;
    double hmax = 905.0  / JG_TO_BTULBM;
    double rtol = 1.0e-6;

    double a = hmin;
    double b = hmax;

    double fa = Temperature(a, p) - T;
    double fb = Temperature(b, p) - T;
    Check(fa * fb < 0.0);

    double h;

    // Iterate until relative tolerance meets stopping criterion
    while (2.0 * std::abs(b - a)/(a + b) > rtol)
    {
        // Evaluate at midpoint
        h = 0.5 * (a + b);
        double fh = Temperature(h, p) - T;

        // If midpoint and left endpoints are of opposite sign,
        // the midpoint becomes new right
        if (fh * fa < 0)
        {
            b = h;
        }
        // Otherwise, midpoint becomes new left
        else
        {
            a = h;
            fa = fh;
        }
    }

    return h;
}

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
// end of Water_Properties.cpp
//---------------------------------------------------------------------------//
