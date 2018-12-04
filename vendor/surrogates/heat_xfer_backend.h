#ifndef MAGNOLIA_HEAT_XFER_BACKEND_H
#define MAGNOLIA_HEAT_XFER_BACKEND_H

void solve_steady_nonlin(double *source, double T_co, double *r_grid_fuel,
  double *r_grid_clad, int n_fuel_rings, int n_clad_rings,
  double tol, double *T);

// From IAPWS
double nu1(double p, double T);

#endif // MAGNOLIA_HEAT_XFER_BACKEND_H
