#include "smrt/Two_Group_Cross_Sections.h"

#include "Nemesis/harness/DBC.hh"

namespace stream
{
//---------------------------------------------------------------------------//
// Get interpolated cross section data
//---------------------------------------------------------------------------//
auto Two_Group_Cross_Sections::get_data(
        Assembly_Model::PIN_TYPE type,
        double                   T,
        double                   rho) -> XS_Data
{
    Require(T > 273);
    Require(T < 3000);
    Require(rho > 0);
    Require(rho < 1);

    constexpr double T_base = 1000;   // degrees K
    constexpr double rho_base = 0.75; // g/cc

    XS_Data data;

    if (type == Assembly_Model::FUEL)
    {
        // XS Data (mean value, temperature slope, density slope
        constexpr double D0[3]       = {0.6109,  -3.45e-7,   -0.5895};
        constexpr double D1[3]       = {0.2164,  -6.795e-07, -0.2568};
        constexpr double siga0[3]    = {0.01607,  5.689e-07,  0.006413};
        constexpr double siga1[3]    = {0.3116,   4.513e-06,  0.03421};
        constexpr double sigs01[3]   = {0.01172, -3.738e-07,  0.02188};
        constexpr double nu_sigf0[3] = {0.01191, -4.358e-08,  0.002904};
        constexpr double nu_sigf1[3] = {0.5007,   5.105e-06,  0.0469};

        auto interp = [T_base, rho_base, T, rho](const double coeffs[2])
        {
            return coeffs[0] +
                   coeffs[1] * (T - T_base) +
                   coeffs[2] * (rho - rho_base);
        };

        // Don't correct D and nu-fission
        data.diffusion[0]  = interp(D0);
        data.diffusion[1]  = interp(D1);
        data.absorption[0] = interp(siga0);
        data.absorption[1] = interp(siga1);
        data.scatter       = interp(sigs01);
        data.nu_fission[0] = nu_sigf0[0];
        data.nu_fission[1] = nu_sigf1[0];
    }
    else if (type == Assembly_Model::GUIDE)
    {
        // XS Data (mean value, temperature slope, density slope
        constexpr double D0[2]       = {0.4916,    -0.7582};
        constexpr double D1[2]       = {0.1718,    -0.2746};
        constexpr double siga0[2]    = {0.0005049,  0.0007154};
        constexpr double siga1[2]    = {0.01996,    0.02726};
        constexpr double sigs01[2]   = {0.03082,    0.04422};

        auto interp = [rho_base, rho](const double coeffs[2])
        {
            return coeffs[0] +
                   coeffs[1] * (rho - rho_base);
        };

        // Don't correct D and nu-fission
        data.diffusion[0]  = D0[0];
        data.diffusion[1]  = D1[0];
        data.absorption[0] = interp(siga0);
        data.absorption[1] = interp(siga1);
        data.scatter       = interp(sigs01);
        data.nu_fission[0] = 0.0;
        data.nu_fission[1] = 0.0;
    }
    return data;
}

//---------------------------------------------------------------------------//
} // end namespace stream

//---------------------------------------------------------------------------//
// end of Two_Group_Cross_Sections.cpp
//---------------------------------------------------------------------------//
