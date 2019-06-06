#include "smrt/Two_Group_Diffusion.h"

// Trilinos stuff
#include "Teuchos_RCP.hpp"
#ifdef COMM_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#else
#include "Teuchos_DefaultSerialComm.hpp"
#endif
#include "BelosSolverFactory.hpp"
#include "BelosTpetraAdapter.hpp"

// vendored includes
#include <gsl/gsl>

namespace enrico {
//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
Two_Group_Diffusion::Two_Group_Diffusion(SP_Assembly assembly,
                                         RCP_PL params,
                                         const Vec_Dbl& dz)
  : d_assembly(assembly)
  , d_dz(dz)
{
  const auto& x_edges = d_assembly->x_edges();
  const auto& y_edges = d_assembly->y_edges();

  d_Nx = x_edges.size() - 1;
  d_Ny = y_edges.size() - 1;
  d_Nz = d_dz.size();
  d_num_cells = d_Nx * d_Ny * d_Nz;

  d_dx.resize(d_Nx);
  for (int ix = 0; ix < d_Nx; ++ix)
    d_dx[ix] = x_edges[ix + 1] - x_edges[ix];

  d_dy.resize(d_Ny);
  for (int iy = 0; iy < d_Ny; ++iy)
    d_dy[iy] = y_edges[iy + 1] - y_edges[iy];

  // Solve parameters
  d_tol = 1e-4;
  d_max_iters = 1000;

  // Process boundary conditions
  d_bcs.resize(6, VACUUM);
  using Array_Str = Teuchos::Array<std::string>;
  if (params->isType<Array_Str>("boundary")) {
    auto bc_str = params->get<Array_Str>("boundary");
    Expects(bc_str.size() == 6);
    for (int bdry = 0; bdry < 6; ++bdry) {
      Expects(bc_str[bdry] == "vacuum" || bc_str[bdry] == "reflect");
      if (bc_str[bdry] == "reflect")
        d_bcs[bdry] = REFLECT;
      else
        d_bcs[bdry] = VACUUM;
    }
  }

  // Make the communicator
#ifdef COMM_MPI
  Teuchos::RCP<const Teuchos::Comm<int>> comm =
    Teuchos::rcp(new Teuchos::MpiComm<int>(nemesis::communicator));
#else
  Teuchos::RCP<const Teuchos::Comm<int>> comm =
    Teuchos::rcp(new Teuchos::SerialComm<int>());
#endif

  // Build map (no decomposition)
  d_map = Teuchos::rcp(new MAP(d_num_cells, d_num_cells, 0, comm));

  // Build matrix graph
  d_graph = Teuchos::rcp(new GRAPH(d_map, 7));
  Vec_Int row_inds;
  row_inds.reserve(7);
  for (int iz = 0; iz < d_Nz; ++iz) {
    for (int iy = 0; iy < d_Ny; ++iy) {
      for (int ix = 0; ix < d_Nx; ++ix) {
        int cell = cellid(ix, iy, iz);

        row_inds.clear();

        // Diagonal -- always present
        row_inds.push_back(cell);

        // Low Z
        if (iz > 0)
          row_inds.push_back(cellid(ix, iy, iz - 1));

        // Low Y
        if (iy > 0)
          row_inds.push_back(cellid(ix, iy - 1, iz));

        // Low X
        if (ix > 0)
          row_inds.push_back(cellid(ix - 1, iy, iz));

        // High X
        if (ix < d_Nx - 1)
          row_inds.push_back(cellid(ix + 1, iy, iz));

        // High y
        if (iy < d_Ny - 1)
          row_inds.push_back(cellid(ix, iy + 1, iz));

        // High z
        if (iz < d_Nz - 1)
          row_inds.push_back(cellid(ix, iy, iz + 1));

        d_graph->insertGlobalIndices(cell, Teuchos::arrayViewFromVector(row_inds));
      }
    }
  }
  d_graph->fillComplete();

  // Create vectors and matrices
  d_A_fast = Teuchos::rcp(new MATRIX(d_graph));
  d_A_therm = Teuchos::rcp(new MATRIX(d_graph));
  d_x_fast = Teuchos::rcp(new VECTOR(d_map));
  d_x_therm = Teuchos::rcp(new VECTOR(d_map));
  d_b_fast = Teuchos::rcp(new VECTOR(d_map));
  d_b_therm = Teuchos::rcp(new VECTOR(d_map));
}

//---------------------------------------------------------------------------//
// Solve problem given fuel temperatures and moderator densities
//---------------------------------------------------------------------------//
void Two_Group_Diffusion::solve(const Vec_Dbl& temperatures,
                                const Vec_Dbl& densities,
                                Vec_Dbl& powers)
{
  Expects(temperatures.size() == d_num_cells);
  Expects(densities.size() == d_num_cells);
  Expects(powers.size() == d_num_cells);

  // Compute cross section data at given conditions
  std::vector<XS_Data> xs_data(d_num_cells);
  for (int iz = 0; iz < d_Nz; ++iz) {
    for (int iy = 0; iy < d_Ny; ++iy) {
      for (int ix = 0; ix < d_Nx; ++ix) {
        int cell = cellid(ix, iy, iz);
        int pin = cellid(ix, iy, 0);

        xs_data[cell] = d_xs.get_data(
          d_assembly->pin_type(ix, iy), temperatures[cell], densities[cell]);
      }
    }
  }

  build_matrices(xs_data);

  // Build Belos solver
  auto belos_pl = Teuchos::rcp(new Teuchos::ParameterList());
  belos_pl->set("Convergence Tolerance", 0.1 * d_tol);
  belos_pl->set("Maximum Iterations", 100);
  belos_pl->set("Verbosity", Belos::Warnings);
  belos_pl->set("Implicit Residual Scaling", "Norm of RHS");
  belos_pl->set("Explicit Residual Scaling", "Norm of RHS");
  Belos::SolverFactory<ST, MV, OP> factory;
  std::string solver_type = "CG";
  auto solver = factory.create(solver_type, belos_pl);

  // Create linear problems
  Teuchos::RCP<Belos::LinearProblem<ST, MV, OP>> fast_problem(
    new Belos::LinearProblem<ST, MV, OP>());
  fast_problem->setOperator(d_A_fast);
  fast_problem->setLHS(d_x_fast);
  fast_problem->setRHS(d_b_fast);

  Teuchos::RCP<Belos::LinearProblem<ST, MV, OP>> therm_problem(
    new Belos::LinearProblem<ST, MV, OP>());
  therm_problem->setOperator(d_A_therm);
  therm_problem->setLHS(d_x_therm);
  therm_problem->setRHS(d_b_therm);

  auto x_fast_data = d_x_fast->getDataNonConst();
  auto x_therm_data = d_x_therm->getDataNonConst();
  auto b_fast_data = d_b_fast->getDataNonConst();
  auto b_therm_data = d_b_therm->getDataNonConst();

  Teuchos::ArrayRCP<double> fission(d_num_cells);
  Teuchos::ArrayRCP<double> fission_diff(d_num_cells);
  std::copy(powers.begin(), powers.end(), fission().begin());

  // Set initial guess for fluxes
  std::fill(x_fast_data.begin(), x_fast_data.end(), 0.0);
  std::fill(x_therm_data.begin(), x_therm_data.end(), 0.0);

  // Normalize, if starting vector is zero then set to constant
  double old_fisn_nrm;
  double fisn_nrm = vec_norm(fission);
  if (fisn_nrm == 0.0) {
    std::fill(fission.begin(), fission.end(), 1.0);
    fisn_nrm = std::sqrt(static_cast<double>(d_num_cells));
  }
  Expects(fisn_nrm > 0.0);
  scale_vec(fission, 1.0 / fisn_nrm);

  // Start solve
  double keff = 1.0;
  double rel_err;
  int it;
  for (it = 0; it < d_max_iters; ++it) {
    old_fisn_nrm = fisn_nrm;
    std::copy(fission.begin(), fission.end(), b_fast_data.begin());

    // Normalize fission source
    scale_vec(b_fast_data, 1.0 / keff);

    // Solve fast group
    fast_problem->setProblem();
    solver->setProblem(fast_problem);
    Belos::ReturnType result = solver->solve();

    // Compute scatter source for thermal group
    for (int cell = 0; cell < d_num_cells; ++cell) {
      b_therm_data[cell] = xs_data[cell].scatter * x_fast_data[cell];
    }

    // Solve thermal group
    therm_problem->setProblem();
    solver->setProblem(therm_problem);
    result = solver->solve();

    // Compute new fission source
    for (int cell = 0; cell < d_num_cells; ++cell) {
      fission[cell] = xs_data[cell].nu_fission[0] * x_fast_data[cell] +
                      xs_data[cell].nu_fission[1] * x_therm_data[cell];
      fission_diff[cell] = fission[cell] - keff * b_fast_data[cell];
    }
    fisn_nrm = vec_norm(fission);
    Expects(fisn_nrm > 0.0);
    rel_err = vec_norm(fission_diff) / fisn_nrm;

    // Compute new keff
    keff *= (fisn_nrm / old_fisn_nrm);
    Expects(keff > 0.0);
    Expects(keff < 2.0);

    // Check for convergence
    if (rel_err < d_tol) {
      std::cout << "Solver converged to keff = " << keff << " after " << it
                << " iterations" << std::endl;
      break;
    }
  }

  if (it == d_max_iters) {
    std::cout << "Solver failed to converge after " << d_max_iters
              << " iterations, final relative residual = " << rel_err << std::endl;
  }

  // Volume-weight and normalize fission source
  for (int iz = 0; iz < d_Nz; ++iz) {
    for (int iy = 0; iy < d_Ny; ++iy) {
      for (int ix = 0; ix < d_Nx; ++ix) {
        int cell = cellid(ix, iy, iz);
        double vol = d_dx[ix] * d_dy[iy] * d_dz[iz];
        fission[cell] *= vol;
      }
    }
  }

  fisn_nrm = vec_norm(fission);
  scale_vec(fission, 1.0 / fisn_nrm);
  std::copy(fission.begin(), fission.end(), powers.begin());
}

//---------------------------------------------------------------------------//
// Build matrices for diffusion equation
//---------------------------------------------------------------------------//
void Two_Group_Diffusion::build_matrices(const std::vector<XS_Data>& xs_data)
{
  Expects(xs_data.size() == d_num_cells);

  // Build matrices
  Vec_Int row_inds;
  Vec_Dbl row_vals_fast, row_vals_therm;
  row_inds.reserve(7);
  row_vals_fast.reserve(7);
  row_vals_therm.reserve(7);
  d_A_fast->resumeFill();
  d_A_therm->resumeFill();
  for (int iz = 0; iz < d_Nz; ++iz) {
    for (int iy = 0; iy < d_Ny; ++iy) {
      for (int ix = 0; ix < d_Nx; ++ix) {
        int cell = cellid(ix, iy, iz);
        int pin = cellid(ix, iy, 0);

        const auto& cell_data = xs_data[cell];

        // Clear matrix row data
        row_inds.clear();
        row_vals_therm.clear();
        row_vals_fast.clear();

        // Diagonal -- always present
        row_inds.push_back(cell);
        row_vals_fast.push_back(0.0);
        row_vals_therm.push_back(0.0);

        // Removal terms
        row_vals_fast[0] = cell_data.absorption[0] + cell_data.scatter;
        row_vals_therm[0] = cell_data.absorption[1];

        // Low Z
        if (iz > 0) {
          int nbr_cell = cellid(ix, iy, iz - 1);
          row_inds.push_back(nbr_cell);

          const auto& nbr_data = xs_data[nbr_cell];

          // Get edge diffusion coefficients
          double D0 = 0.5 * (cell_data.diffusion[0] + nbr_data.diffusion[0]);
          double D1 = 0.5 * (cell_data.diffusion[1] + nbr_data.diffusion[1]);

          // Average delta
          double delta = 0.5 * (d_dz[iz] + d_dz[iz - 1]);

          double val0 = D0 / (delta * delta);
          double val1 = D1 / (delta * delta);
          row_vals_fast.push_back(-val0);
          row_vals_fast[0] += val0;
          row_vals_therm.push_back(-val1);
          row_vals_therm[0] += val1;
        } else if (d_bcs[LO_Z] == VACUUM) {
          row_vals_fast[0] += 2.0 * cell_data.diffusion[0] / (d_dz[iz] * d_dz[iz]);
          row_vals_therm[0] += 2.0 * cell_data.diffusion[1] / (d_dz[iz] * d_dz[iz]);
        }

        // Low Y
        if (iy > 0) {
          int nbr_cell = cellid(ix, iy - 1, iz);
          row_inds.push_back(nbr_cell);

          const auto& nbr_data = xs_data[nbr_cell];

          // Get edge diffusion coefficients
          double D0 = 0.5 * (cell_data.diffusion[0] + nbr_data.diffusion[0]);
          double D1 = 0.5 * (cell_data.diffusion[1] + nbr_data.diffusion[1]);

          // Average delta
          double delta = 0.5 * (d_dy[iy] + d_dy[iy - 1]);

          double val0 = D0 / (delta * delta);
          double val1 = D1 / (delta * delta);
          row_vals_fast.push_back(-val0);
          row_vals_fast[0] += val0;
          row_vals_therm.push_back(-val1);
          row_vals_therm[0] += val1;
        } else if (d_bcs[LO_Y] == VACUUM) {
          row_vals_fast[0] += 2.0 * cell_data.diffusion[0] / (d_dy[iy] * d_dy[iy]);
          row_vals_therm[0] += 2.0 * cell_data.diffusion[1] / (d_dy[iy] * d_dy[iy]);
        }

        // Low X
        if (ix > 0) {
          int nbr_cell = cellid(ix - 1, iy, iz);
          row_inds.push_back(nbr_cell);

          const auto& nbr_data = xs_data[nbr_cell];

          // Get edge diffusion coefficients
          double D0 = 0.5 * (cell_data.diffusion[0] + nbr_data.diffusion[0]);
          double D1 = 0.5 * (cell_data.diffusion[1] + nbr_data.diffusion[1]);

          // Average delta
          double delta = 0.5 * (d_dx[ix] + d_dx[ix - 1]);

          double val0 = D0 / (delta * delta);
          double val1 = D1 / (delta * delta);
          row_vals_fast.push_back(-val0);
          row_vals_fast[0] += val0;
          row_vals_therm.push_back(-val1);
          row_vals_therm[0] += val1;
        } else if (d_bcs[LO_X] == VACUUM) {
          row_vals_fast[0] += 2.0 * cell_data.diffusion[0] / (d_dx[ix] * d_dx[ix]);
          row_vals_therm[0] += 2.0 * cell_data.diffusion[1] / (d_dx[ix] * d_dx[ix]);
        }

        // High X
        if (ix < d_Nx - 1) {
          int nbr_cell = cellid(ix + 1, iy, iz);
          row_inds.push_back(nbr_cell);

          const auto& nbr_data = xs_data[nbr_cell];

          // Get edge diffusion coefficients
          double D0 = 0.5 * (cell_data.diffusion[0] + nbr_data.diffusion[0]);
          double D1 = 0.5 * (cell_data.diffusion[1] + nbr_data.diffusion[1]);

          // Average delta
          double delta = 0.5 * (d_dx[ix] + d_dx[ix + 1]);

          double val0 = D0 / (delta * delta);
          double val1 = D1 / (delta * delta);
          row_vals_fast.push_back(-val0);
          row_vals_fast[0] += val0;
          row_vals_therm.push_back(-val1);
          row_vals_therm[0] += val1;
        } else if (d_bcs[HI_X] == VACUUM) {
          row_vals_fast[0] += 2.0 * cell_data.diffusion[0] / (d_dx[ix] * d_dx[ix]);
          row_vals_therm[0] += 2.0 * cell_data.diffusion[1] / (d_dx[ix] * d_dx[ix]);
        }

        // High y
        if (iy < d_Ny - 1) {
          int nbr_cell = cellid(ix, iy + 1, iz);
          row_inds.push_back(nbr_cell);

          const auto& nbr_data = xs_data[nbr_cell];

          // Get edge diffusion coefficients
          double D0 = 0.5 * (cell_data.diffusion[0] + nbr_data.diffusion[0]);
          double D1 = 0.5 * (cell_data.diffusion[1] + nbr_data.diffusion[1]);

          // Average delta
          double delta = 0.5 * (d_dy[iy] + d_dy[iy + 1]);

          double val0 = D0 / (delta * delta);
          double val1 = D1 / (delta * delta);
          row_vals_fast.push_back(-val0);
          row_vals_fast[0] += val0;
          row_vals_therm.push_back(-val1);
          row_vals_therm[0] += val1;
        } else if (d_bcs[HI_Y] == VACUUM) {
          row_vals_fast[0] += 2.0 * cell_data.diffusion[0] / (d_dy[iy] * d_dy[iy]);
          row_vals_therm[0] += 2.0 * cell_data.diffusion[1] / (d_dy[iy] * d_dy[iy]);
        }

        // High z
        if (iz < d_Nz - 1) {
          int nbr_cell = cellid(ix, iy, iz + 1);
          row_inds.push_back(nbr_cell);

          const auto& nbr_data = xs_data[nbr_cell];

          // Get edge diffusion coefficients
          double D0 = 0.5 * (cell_data.diffusion[0] + nbr_data.diffusion[0]);
          double D1 = 0.5 * (cell_data.diffusion[1] + nbr_data.diffusion[1]);

          // Average delta
          double delta = 0.5 * (d_dz[iz] + d_dz[iz + 1]);

          double val0 = D0 / (delta * delta);
          double val1 = D1 / (delta * delta);
          row_vals_fast.push_back(-val0);
          row_vals_fast[0] += val0;
          row_vals_therm.push_back(-val1);
          row_vals_therm[0] += val1;
        } else if (d_bcs[HI_Z] == VACUUM) {
          row_vals_fast[0] += 2.0 * cell_data.diffusion[0] / (d_dz[iz] * d_dz[iz]);
          row_vals_therm[0] += 2.0 * cell_data.diffusion[1] / (d_dz[iz] * d_dz[iz]);
        }

        // Add values to matrices
        d_A_fast->replaceGlobalValues(cell,
                                      Teuchos::arrayViewFromVector(row_inds),
                                      Teuchos::arrayViewFromVector(row_vals_fast));
        d_A_therm->replaceGlobalValues(cell,
                                       Teuchos::arrayViewFromVector(row_inds),
                                       Teuchos::arrayViewFromVector(row_vals_therm));

      } // ix
    }   // iy
  }     // iz

  // Complete matrices
  d_A_fast->fillComplete();
  d_A_therm->fillComplete();
}

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
// end of Two_Group_Diffusion.cpp
//---------------------------------------------------------------------------//
