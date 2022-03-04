#include "enrico/boron_driver.h"

#include <gsl/gsl>  // For Expects

namespace enrico {

BoronDriver::BoronDriver(MPI_Comm comm, pugi::xml_node node)
  : Driver(comm)
{
  MPI_Barrier(MPI_COMM_WORLD);

  // Get the target k_eff to search for
  if (node.child("target_keff")) {
    target_k_eff_ = node.child("target_keff").text().as_double();
  } else {
    target_k_eff_ = 1.;
  }
  Expects(target_k_eff_ > 0.);

  // Get the B10 isotopic abundance to use, if available
  if (node.child("B10_enrichment")) {
    B10_iso_abund_ = node.child("B10_enrichment").text().as_double();
  }
  Expects(B10_iso_abund_ > 0.);
  Expects(B10_iso_abund_ <= 1.);

  // Get the convergence epsilon to use when checking for convergence
  if (node.child("boron_epsilon")) {
    epsilon_ = nodechild("boron_epsilon").text().as_double();
  }
  Expects(epsilon_ > 0.);
  Expects(epsilon_ < 1.);
}

double BoronDriver::solve_ppm(bool first_pass, double k_eff, double k_eff_prev)
{
  // The estimation is performed using the Newton-Raphson method
  double del_ppm_del_k;
  double new_ppm;
  if (first_pass) {
    // Make a reasonable update so we can get two points to compute a derivative
    if (k_eff > target_k_eff_) {
      new_ppm = ppm_prev_ + 1000.;
    } else if (k_eff < target_k_eff_) {
      new_ppm = ppm_prev_ - 1000.;
    } else {
      // Then our most recent calculation was on target, so we keep it
      // This should generally never be achieved by a realistic CE MC case due
      // to the noise, but, it may be reached with simpler MG MC cases or
      // deterministic solvers
      new_ppm = ppm_prev_;
    }

  } else {
    del_ppm_del_k = (ppm_ - ppm_prev_) / (k_eff - k_eff_prev);
    new_ppm = (target_k_eff_ - k_eff) * del_ppm_del_k + ppm_;
  }

  // With noisy results we could have a negative guess for the ppm; this is
  // clearly non-physical. Additionally the first guess value may be negative.
  // For both negative cases a floor of 0 ppm is set
  if (new_ppm < 0.) {
    new_ppm = 0.;
  }

  ppm_prev_ = ppm_;
  ppm_ = new_ppm;

  return ppm_;
}

bool BoronDriver::is_converged() {
  bool converged;
  double norm;

  norm = ppm_ > 0. ? (ppm_ - ppm_prev_) / ppm_: (ppm_ - ppm_prev_);
  converged = norm < epsilon_;

  return converged
}

void BoronDriver::set_fluid_cells(std::vector<CellHandle>& fluid_cell_handles) {
  fluid_cell_handles_ = fluid_cell_handles;
}

void BoronDriver::print_boron()
{
  std::ios_base::fmtflags old_flags(std::cout.flags());
  std::stringstream msg;
  msg << "\tUpdated Boron Concentration [ppm]: " << ppm_;
  comm_.message(msg.str());
  msg.clear();
  msg.str(std::string());
  msg << "\t         B-10 Concentration [ppm]: " << B10_iso_abund_ * ppm_;
  comm_.message(msg.str());
  msg.clear();
  msg.str(std::string());

  std::cout.flags(old_flags);
}
BoronDriver::~BoronDriver()
{

  MPI_Barrier(MPI_COMM_WORLD);
}
} // namespace enrico
