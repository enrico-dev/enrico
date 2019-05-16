#include "iapws/iapws.h"
#include <cmath>

namespace iapws {

//==============================================================================

double gamma1(double pi, double tau)
{
  double out = 0.0;
  for (int i = 0; i < 34; i++) {
    out += _n1f[i] * pow(7.1 - pi, _I1f[i]) * pow(tau - 1.222, _J1f[i]);
  }
  return out;
}

//==============================================================================

double gamma1_pi(double pi, double tau)
{
  double out = 0.0;
  for (int i = 0; i < 34; i++) {
    out -= _n1f[i] * _I1f[i] * pow(7.1 - pi, _I1f[i] - 1) * pow(tau - 1.222, _J1f[i]);
  }
  return out;
}

//==============================================================================

double gamma1_tau(double pi, double tau)
{
  double out = 0.0;
  for (int i = 0; i < 34; i++) {
    out += _n1f[i] * pow(7.1 - pi, _I1f[i]) * _J1f[i] * pow(tau - 1.222, _J1f[i] - 1);
  }
  return out;
}

//==============================================================================

double gamma1_tau_tau(double pi, double tau)
{
  double out = 0.0;
  for (int i = 0; i < 34; i++) {
    out += _n1f[i] * pow(7.1 - pi, _I1f[i]) * _J1f[i] * (_J1f[i] - 1) *
           pow(tau - 1.222, _J1f[i] - 2);
  }
  return out;
}

//==============================================================================
// Region 1 thermodynamic properties.
//==============================================================================

double nu1(double p, double T)
{
  double pi = p / 16.53;   // Dimensionless pressure; p must be in MPa.
  double tau = 1386.0 / T; // Dimensionless temperature; T must be in K.
  return pi * gamma1_pi(pi, tau) * R * T / (p * 1e3);
}

//==============================================================================

double u1(double p, double T)
{
  double pi = p / 16.53;   // Dimensionless pressure; p must be in MPa.
  double tau = 1386.0 / T; // Dimensionless temperature; T must be in K.
  return (tau * gamma1_tau(pi, tau) - pi * gamma1_pi(pi, tau)) * R * T;
}

//==============================================================================

double s1(double p, double T)
{
  double pi = p / 16.53;   // Dimensionless pressure; p must be in MPa.
  double tau = 1386.0 / T; // Dimensionless temperature; T must be in K.
  return (tau * gamma1_tau(pi, tau) - gamma1(pi, tau)) * R;
}

//==============================================================================

double h1(double p, double T)
{
  double pi = p / 16.53;   // Dimensionless pressure; p must be in MPa.
  double tau = 1386.0 / T; // Dimensionless temperature; T must be in K.
  return tau * gamma1_tau(pi, tau) * R * T;
}

//==============================================================================

double cp1(double p, double T)
{
  double pi = p / 16.53;   // Dimensionless pressure; p must be in MPa.
  double tau = 1386.0 / T; // Dimensionless temperature; T must be in K.
  return -tau * tau * gamma1_tau_tau(pi, tau) * R;
}

//==============================================================================
// Region 1 backward equations.
//==============================================================================

double T_from_p_h(double p, double h)
{
  double pi = p;           // p must be in MPa
  double eta = h / 2500.0; // h must be in kJ / kg
  double theta = 0.0;
  for (int i = 0; i < 20; i++) {
    theta += _n1bh[i] * pow(pi, _I1bh[i]) * pow(eta + 1, _J1bh[i]);
  }
  return theta;
}

//==============================================================================
// Region 2 forward equations.
//==============================================================================

double gamma2ig(double pi, double tau)
{
  double out = log(pi);
  for (int i = 0; i < 9; i++) {
    out += _n2igf[i] * pow(tau, _J2igf[i]);
  }
  return out;
}

//==============================================================================

double gamma2r(double pi, double tau)
{
  double out = 0.0;
  for (int i = 0; i < 43; i++) {
    out += _n2rf[i] * pow(pi, _I2rf[i]) * pow(tau - 0.5, _J2rf[i]);
  }
  return out;
}

//==============================================================================

double gamma2ig_pi(double pi, double tau)
{
  return 1.0 / pi;
}

//==============================================================================

double gamma2r_pi(double pi, double tau)
{
  double out = 0.0;
  for (int i = 0; i < 43; i++) {
    out += _n2rf[i] * _I2rf[i] * pow(pi, _I2rf[i] - 1) * pow(tau - 0.5, _J2rf[i]);
  }
  return out;
}

//==============================================================================

double gamma2ig_tau(double pi, double tau)
{
  double out = 0.0;
  for (int i = 0; i < 9; i++) {
    out += _n2igf[i] * _J2igf[i] * pow(tau, _J2igf[i] - 1);
  }
  return out;
}

//==============================================================================

double gamma2r_tau(double pi, double tau)
{
  double out = 0.0;
  for (int i = 0; i < 43; i++) {
    out += _n2rf[i] * pow(pi, _I2rf[i]) * _J2rf[i] * pow(tau - 0.5, _J2rf[i] - 1);
  }
  return out;
}

//==============================================================================
// Region 2 thermodynamic properties.
//==============================================================================

double nu2(double p, double T)
{
  double pi = p;          // Dimensionless pressure; p must be in MPa.
  double tau = 540.0 / T; // Dimensionless temperature; T must be in K.
  return pi * (gamma2ig_pi(pi, tau) + gamma2r_pi(pi, tau)) * R * T / (p * 1e3);
}

//==============================================================================

double h2(double p, double T)
{
  double pi = p;          // Dimensionless pressure; p must be in MPa.
  double tau = 540.0 / T; // Dimensionless temperature; T must be in K.
  return tau * (gamma2ig_tau(pi, tau) + gamma2r_tau(pi, tau)) * R * T;
}

//==============================================================================
// Region 4 (saturation line) equations.
//==============================================================================

double sat_press(double T)
{
  double theta = T + _n4[8] / (T - _n4[9]);
  double A = pow(theta, 2) + _n4[0] * theta + _n4[1];
  double B = _n4[2] * pow(theta, 2) + _n4[3] * theta + _n4[4];
  double C = _n4[5] * pow(theta, 2) + _n4[6] * theta + _n4[7];
  return pow(2.0 * C / (-B + sqrt(pow(B, 2) - 4 * A * C)), 4);
}

//==============================================================================

double sat_temp(double p)
{
  double beta = pow(p, 0.25);
  double E = pow(beta, 2) + _n4[2] * beta + _n4[5];
  double F = _n4[0] * pow(beta, 2) + _n4[3] * beta + _n4[6];
  double G = _n4[1] * pow(beta, 2) + _n4[4] * beta + _n4[7];
  double D = 2.0 * G / (-F - sqrt(pow(F, 2) - 4 * E * G));
  return 0.5 * (_n4[9] + D - sqrt(pow(_n4[9] + D, 2) - 4.0 * (_n4[8] + _n4[9] * D)));
}

//==============================================================================
// Viscosity equations.
//==============================================================================

double mu_0(double T_bar)
{
  double denominator = 0.0;
  for (int i = 0; i < 4; i++) {
    denominator += _H_i_vals[i] / pow(T_bar, i);
  }
  return 100.0 * sqrt(T_bar) / denominator;
}

//==============================================================================

double mu_1(double T_bar, double rho_bar)
{
  double out = 0.0;
  for (int i = 0; i < 6; i++) {
    double coeff = 0.0;
    for (int j = 0; j < 7; j++) {
      coeff += _H_ij_vals[7 * i + j] * pow(rho_bar - 1, j);
    }
    out += coeff * pow(1.0 / T_bar - 1, i);
  }
  return exp(rho_bar * out);
}

//==============================================================================

double mu(double T, double rho)
{
  double T_bar = T / 647.096;
  double rho_bar = rho / 322.0;
  return mu_0(T_bar) * mu_1(T_bar, rho_bar) * 1e-6;
}

//==============================================================================
// Thermal conductivity equations.
//==============================================================================

double k_0(double T_bar)
{
  double denominator = 0.0;
  for (int k = 0; k < 5; k++) {
    denominator += _L_k_vals[k] / pow(T_bar, k);
  }
  return sqrt(T_bar) / denominator;
}

//==============================================================================

double k_1(double T_bar, double rho_bar)
{
  double out = 0.0;
  for (int i = 0; i < 5; i++) {
    double coeff = 0.0;
    for (int j = 0; j < 6; j++) {
      coeff += _L_ij_vals[6 * i + j] * pow(rho_bar - 1, j);
    }
    out += coeff * pow(1.0 / T_bar - 1, i);
  }
  return exp(rho_bar * out);
}

//==============================================================================

double k(double T, double rho)
{
  double T_bar = T / 647.096;
  double rho_bar = rho / 322.0;
  return k_0(T_bar) * k_1(T_bar, rho_bar) * 1e-3;
}

//==============================================================================
// Surface tension.
//==============================================================================

double surf_tens(double T)
{
  double tau = 1.0 - T / 647.096;
  return 0.2358 * pow(tau, 1.256) * (1.0 - 0.625 * tau);
}

}
