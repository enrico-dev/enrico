#include <cmath>
#include <iostream>
#include <stdlib.h>

#include "iapws.cpp"

const double PI = 3.1415926535898;

// Declare TH parameters.
const double K_CLAD = 13.0;  // T&K pg 394
//const double H_GAP = 4300.0;  // T&K pg 428
const double RHO_FUEL = 10.97e3;  // kg/m^3 T&K pg 392
const double RHO_CLAD = 6.5e3;  // kg/m^3 T&K pg 394
const double CP_CLAD = 330;  // J/kg/K T&K pg 394

//==============================================================================
// k_fuel
//==============================================================================

double
k_fuel(double T)
{
  // T&K eq 8.22c, pg 402.
  double A = 4.52;  // cm K / W
  double B = 2.46e-2;  // cm / W
  double E = 3.5e7;  // W K / cm
  double F = 16361;  // K
  double k = 1 / (A + B*T) + E / (T * T) * std::exp(-F / T);  // W / cm / K
  return k * 100;  // Converted to W / m / K
}

//==============================================================================
// cp_fuel
//==============================================================================

double
cp_fuel(double T)
{
  // T&K eq 8.27, pg 375.
  double t = T / 1000;
  double c2 = 193.238;
  double c3 = 162.8647;
  double c4 = -104.0014;
  double c5 = 29.2056;
  double c6 = -1.9507;
  double c7 = 2.6441;

  double t2 = t * t;
  double t3 = t2 * t;
  double t4 = t3 * t;

  double cp = c2 + 2*c3*t + 3*c4*t2 + 4*c5*t3 + 5*c6*t4 - c7/t2;
  return cp;
}

//==============================================================================
// emiss_fuel
//==============================================================================

double
emiss_fuel(double T)
{
  // MATPRO-11 section 3.
  if (T >= 1000.0) {
    if (T <= 2050.0) {
      return 1.311 - 4.404e-4 * T;
    } else {
      return 0.4083;
    }
  } else {
    return 0.8707;
  }
}

//==============================================================================
// gap_rad_htc
//==============================================================================

double
gap_rad_htc(double T_fo, double T_ci, double r_fo, double r_ci)
{
  double stefan_boltzmann = 5.6697e-8;  // [W / m^2 / K^4]

  double e_fuel = emiss_fuel(T_fo);
  double e_clad = 0.5;  // TODO: verify this

  double F = 1.0 / (1.0/e_fuel + (r_fo / r_ci) * (1.0 / e_clad - 1));

  return stefan_boltzmann * F * (T_fo*T_fo + T_ci*T_ci) * (T_fo + T_ci);
}

//==============================================================================
// gap_gas_htc
//==============================================================================

double
gap_gas_htc(double T, double gap_thick)
{
  double k = 3.366e-3 * std::pow(T, 0.668);  // Assuming only He gas, from MATPRO-11
  double eff_thick = 10.0e-6 + gap_thick;  // T&K pg 419
  return k / eff_thick;
}

//==============================================================================
// gap_htc
//==============================================================================

double
gap_htc(double T_fo, double T_ci, double r_fo, double r_ci)
{
  double T_avg = (T_fo + T_ci) / 2.0;
  double thickness = r_ci - r_fo;
  return gap_rad_htc(T_fo, T_ci, r_fo, r_ci) + gap_gas_htc(T_avg, thickness);
}

//==============================================================================
// chen_htc_nb_const
//==============================================================================

double
chen_htc_nb_const(double Re, double P, double T_b)
{
  double T_sat = sat_temp(P);
  double S = 1.0 / (1.0 + 2.53e-6 * std::pow(Re, 1.17));
  double rho_f = 1.0 / nu1(P, T_sat);
  double k_f = k(T_sat, rho_f);
  double cp_f = cp1(P, T_sat) * 1e3;
  double mu_f = mu(T_sat, rho_f);
  double h_fg = (h2(P, T_sat) - h1(P, T_sat)) *1e3;
  double rho_g = 1.0 / nu2(P, T_sat);
  double sigma = surf_tens(T_b);
  return (S * 0.00122 * std::pow(k_f, 0.79) * std::pow(cp_f, 0.45) * std::pow(rho_f, 0.49)
          / (std::pow(sigma, 0.5) * std::pow(mu_f, 0.29) * std::pow(h_fg, 0.25)
             * std::pow(rho_g, 0.24)));
}

//==============================================================================
//==============================================================================

std::pair<double, double>
find_chen_nb_htc(double heat_flux, double htc_conv, double chen_const,
                 double T_b, double T_sat, double P, double tol=1e-4)
{
  // Assume the cladding temperature is close to saturation.  Raise the
  // temperature in 1 degree increments until we overshoot the answer.
  double htc_nb;
  double T_co_lhs = INFINITY;
  double T_co_rhs = T_sat + 0.1 - 1.0;
  while (T_co_lhs > T_co_rhs) {
    T_co_rhs += 1.0;
    double deltaT_sat = T_co_rhs - T_sat;
    double deltaP = (sat_press(T_co_rhs) - P) * 1e6;
    htc_nb = chen_const * std::pow(deltaP, 0.75) * std::pow(deltaT_sat, 0.24);
    T_co_lhs = ((htc_conv * T_b + htc_nb * T_sat + heat_flux)
                / (htc_conv + htc_nb));
  }

  // Use a binary search to find the cladding temperature.
  double T_lo = std::max(T_sat + 0.1, T_co_rhs - 1.0);
  double T_hi = T_co_rhs;
  double last_htc = htc_nb;
  while (true) {
    double T_mid = 0.5 * (T_lo + T_hi);
    double deltaT_sat = T_mid - T_sat;
    double deltaP = (sat_press(T_mid) - P) * 1e6;
    htc_nb = chen_const * std::pow(deltaP, 0.75) * std::pow(deltaT_sat, 0.24);
    double T_co = ((htc_conv * T_b + htc_nb * T_sat + heat_flux)
                   / (htc_conv + htc_nb));
    if (std::abs((htc_nb - last_htc) / last_htc) < tol) {
      return {T_mid, htc_nb};
    } else if (T_co < T_mid) {
      T_hi = T_mid;
    } else {
      T_lo = T_mid;
    }
    last_htc = htc_nb;
  }
}

//==============================================================================
// fill_matrix
//==============================================================================

void
fill_matrix(double *T, double *r_grid_fuel, double *r_grid_clad,
            int n_fuel_rings, int n_clad_rings, double *out)
{
  int n_col = n_fuel_rings + n_clad_rings + 1;

  // Compute the gap htc.
  // TODO: use the temperature on the surface of the fuel and clad rather than
  // the bounding rings.
  double h_gap = gap_htc(T[n_fuel_rings-1], T[n_fuel_rings],
                         r_grid_fuel[n_fuel_rings], r_grid_clad[0]);

  // Declare intermediate variables.
  double r1, r2, r3, r4, r5, r_avg, k_l, k_c, k_r, dr_l, dr_c, dr_r, FD_l, FD_r;

  //============================================================================
  // Fill the matrix terms for the fuel.
  //============================================================================

  // Handle the terms for the inner-most fuel ring.  Note that there is no
  // leakage term on the inward side of this ring.
  r2 = r_grid_fuel[0];
  r3 = r_grid_fuel[1];
  r4 = r_grid_fuel[2];
  r_avg = (r2 + r3) / 2.0;
  k_c = k_fuel(T[0]);
  k_r = k_fuel(T[1]);
  dr_c = r3 - r2;
  dr_r = r4 - r3;
  FD_r = 2.0 * k_c * k_r / (dr_c*k_r + dr_r*k_c) * r3 / r_avg / dr_c;
  out[n_col*0 + 1] = -FD_r;  // Upper diagonal
  out[n_col*1 + 0] = FD_r;   // Main diagonal

  // Iterate over fuel rings.
  for (int i = 1; i < n_fuel_rings - 1; i++)
  {
    // Get the radii of the relevant grid points.
    r1 = r_grid_fuel[i-1];
    r2 = r_grid_fuel[i];
    r3 = r_grid_fuel[i+1];
    r4 = r_grid_fuel[i+2];

    // Compute the average radius of this ring.
    r_avg = (r2 + r3) / 2.0;

    // Get the conductivity and widths of this ring and its neighbors.
    k_l = k_c;
    k_c = k_r;
    k_r = k_fuel(T[i+1]);
    dr_l = r2 - r1;
    dr_c = r3 - r2;
    dr_r = r4 - r3;

    // Compute the finite differencing coefficients and modify the matrix.
    FD_l = 2.0 * k_l * k_c / (dr_l*k_c + dr_c*k_l) * r2 / r_avg / dr_c;
    FD_r = 2.0 * k_c * k_r / (dr_c*k_r + dr_r*k_c) * r3 / r_avg / dr_c;
    out[n_col*0 + i+1] = -FD_r;      // Upper diagonal
    out[n_col*1 + i] = FD_l + FD_r;  // Main diagonal
    out[n_col*2 + i-1] = -FD_l;      // Lower diagonal
  }

  // Handle the terms for the outer-most fuel ring.
  r1 = r_grid_fuel[n_fuel_rings-2];
  r2 = r_grid_fuel[n_fuel_rings-1];
  r3 = r_grid_fuel[n_fuel_rings];
  r4 = r_grid_clad[0];
  r5 = r_grid_clad[1];
  r_avg = (r2 + r3) / 2.0;
  k_l = k_c;
  k_c = k_r;
  k_r = K_CLAD;
  dr_l = r2 - r1;
  dr_c = r3 - r2;
  dr_r = r5 - r4;
  FD_l = 2.0 * k_l * k_c / (dr_l*k_c + dr_c*k_l) * r2 / r_avg / dr_c;
  FD_r = (1.0
          / (r3 / r4 / h_gap + dr_c / 2.0 / k_c + r3 / r4 * dr_r / 2.0 / k_r)
          * r3 / r_avg / dr_c);
  out[n_col*0 + n_fuel_rings] = -FD_r;          // Upper diagonal
  out[n_col*1 + n_fuel_rings-1] = FD_l + FD_r;  // Main diagonal
  out[n_col*2 + n_fuel_rings-2] = -FD_l;        // Lower diagonal

  //============================================================================
  // Fill the matrix terms for the clad.
  //============================================================================

  // Handle the terms for the inner-most clad ring.
  r1 = r_grid_fuel[n_fuel_rings-1];
  r2 = r_grid_fuel[n_fuel_rings];
  r3 = r_grid_clad[0];
  r4 = r_grid_clad[1];
  r5 = r_grid_clad[2];
  r_avg = (r3 + r4) / 2.0;
  k_l = k_c;
  k_c = K_CLAD;
  k_r = K_CLAD;
  dr_l = r2 - r1;
  dr_c = r4 - r3;
  dr_r = r5 - r4;
  FD_l = (1.0
          / (1.0 / h_gap + r3 / r2 * dr_l / 2.0 / k_l + dr_c / 2.0 / k_c)
          * r3 / r_avg / dr_c);
  FD_r = 2.0 * k_c * k_r / (dr_c*k_r + dr_r*k_c) * r4 / r_avg / dr_c;
  out[n_col*0 + n_fuel_rings+1] = -FD_r;      // Upper diagonal
  out[n_col*1 + n_fuel_rings] = FD_l + FD_r;  // Main diagonal
  out[n_col*2 + n_fuel_rings-1] = -FD_l;      // Lower diagonal

  // Iterate over clad rings.
  for (int i = 1; i < n_clad_rings - 1; i++)
  {
    // Get the radii of the relevant grid points.
    r1 = r_grid_clad[i-1];
    r2 = r_grid_clad[i];
    r3 = r_grid_clad[i+1];
    r4 = r_grid_clad[i+2];

    // Compute the average radius of this ring.
    r_avg = (r2 + r3) / 2.0;

    // Get the conductivity and widths of this ring and its neighbors.
    k_l = K_CLAD;
    k_c = K_CLAD;
    k_r = K_CLAD;
    dr_l = r2 - r1;
    dr_c = r3 - r2;
    dr_r = r4 - r3;

    // Compute the finite differencing coefficients and modify the matrix.
    FD_l = 2.0 * k_l * k_c / (dr_l*k_c + dr_c*k_l) * r2 / r_avg / dr_c;
    FD_r = 2.0 * k_c * k_r / (dr_c*k_r + dr_r*k_c) * r3 / r_avg / dr_c;
    out[n_col*0 + n_fuel_rings+i+1] = -FD_r;      // Upper diagonal
    out[n_col*1 + n_fuel_rings+i] = FD_l + FD_r;  // Main diagonal
    out[n_col*2 + n_fuel_rings+i-1] = -FD_l;      // Lower diagonal
  }

  // Set the main diagonal term for the outer-most clad ring.
  r1 = r_grid_clad[n_clad_rings-2];
  r2 = r_grid_clad[n_clad_rings-1];
  r3 = r_grid_clad[n_clad_rings];
  r_avg = (r2 + r3) / 2.0;
  k_l = K_CLAD;
  k_c = K_CLAD;
  dr_l = r2 - r1;
  dr_c = r3 - r2;
  FD_l = 2.0 * k_l * k_c / (dr_l*k_c + dr_c*k_l) * r2 / r_avg / dr_c;
  out[n_col*1 + n_fuel_rings+n_clad_rings-1] = FD_l;   // Main diagonal
  out[n_col*2 + n_fuel_rings+n_clad_rings-2] = -FD_l;  // Lower diagonal

  //============================================================================
  // Add a matrix term for the pseudo-node representing the coolant temperature.
  //============================================================================

  out[n_col*1 + n_fuel_rings+n_clad_rings] = 1;
}

//==============================================================================
// Boundary conditions
//==============================================================================

void
add_dirichlet_bc(double *r_grid_fuel, double *r_grid_clad, int n_fuel_rings,
                 int n_clad_rings, double *A)
{
  int n_col = n_fuel_rings + n_clad_rings + 1;

  // Compute the differencing coefficient.
  double r_mid = (r_grid_clad[n_clad_rings-1] + r_grid_clad[n_clad_rings])
                 / 2.0;
  double k_c = K_CLAD;
  double dr_c = r_grid_clad[n_clad_rings] - r_grid_clad[n_clad_rings-1];
  double FD_r = 2.0 * k_c / dr_c * r_grid_clad[n_clad_rings] / r_mid / dr_c;

  // Set the upper diagonal term.
  A[n_col*0 + n_fuel_rings+n_clad_rings] = -FD_r;

  // Modify the main diagonal term.
  A[n_col*1 + n_fuel_rings+n_clad_rings-1] += FD_r;
}

//==============================================================================
// Heat system matrix solver
//==============================================================================

void
solve_heat_system(double *A, double T_b, double *source, int n_cols, double *T)
{
  int n_rings = n_cols - 1;

  A[n_cols*0 + 1] /= A[n_cols*1 + 0];
  T[0] = source[0] / A[n_cols*1 + 0];
  for (int i = 1; i < n_rings; i++)
  {
    A[n_cols*0 + i+1] /= A[n_cols*1 + i] - A[n_cols*2 + i-1] * A[n_cols*0 + i];
    T[i] = (source[i] - A[n_cols*2 + i-1] * T[i-1])
           / (A[n_cols*1 + i] - A[n_cols*2 + i-1] * A[n_cols*0 + i]);
  }

  T[n_rings-1] -= A[n_cols*0 + n_rings] * T_b;
  for (int i = n_rings-2; i > -1; i--)
  {
    T[i] -= A[n_cols*0 + i+1] * T[i+1];
  }
}

//==============================================================================
//
//==============================================================================

void
solve_steady_nonlin(double *source, double T_co, double *r_grid_fuel,
                    double *r_grid_clad, int n_fuel_rings, int n_clad_rings,
                    double tol, double *T)
{
  // Get the number of columns in the matrix.
  int n_rings = n_fuel_rings + n_clad_rings;
  int n_cols = n_rings + 1;
  double A[3*n_cols];
  double T_last_k_iter[n_rings];

  // Save the current temperature distribution.
  for (int i = 0; i < n_fuel_rings + n_clad_rings; i++) {
    T_last_k_iter[i] = T[i];
  }

  // Iterate on fuel thermal conductivity.
  bool converged = false;
  int iteration = 0;
  while (!converged) {
    iteration++;

    // Solve the linearized system for a new temperature distribution.
    fill_matrix(T, r_grid_fuel, r_grid_clad, n_fuel_rings, n_clad_rings, A);
    add_dirichlet_bc(r_grid_fuel, r_grid_clad, n_fuel_rings, n_clad_rings, A);
    solve_heat_system(A, T_co, source, n_cols, T);

    // Compute the L2 temperature error and check convergence.
    double l2_err = 0.0;
    for (int i = 0; i < n_fuel_rings + n_clad_rings; i++) {
      if (T_last_k_iter[i] != 0) {
        double rel_diff = (T[i] - T_last_k_iter[i]) / T_last_k_iter[i];
        l2_err += rel_diff * rel_diff;
      }
    }
    l2_err = std::sqrt(l2_err);
    converged = l2_err < tol;

    // Save this temperature distribution for convergence checking.
    for (int i = 0; i < n_fuel_rings + n_clad_rings; i++) {
      T_last_k_iter[i] = T[i];
    }
  }
}
