#include "enrico/boron_driver.h"
#include "enrico/const.h"

#include <gsl/gsl>  // For Expects

namespace enrico {

BoronDriver::BoronDriver(MPI_Comm comm, pugi::xml_node node)
  : Driver(comm)
{
  MPI_Barrier(MPI_COMM_WORLD);

  // Get the target k_eff to search for
  if (node.child("target_keff")) {
    target_k_eff_ = node.child("target_keff").text().as_double();
  }
  Expects(target_k_eff_ > 0.);

  if (node.child("target_keff_tolerance")) {
    target_k_eff_tol_ = node.child("target_keff_tolerance").text().as_double();
  }
  Expects(target_k_eff_tol_ > 0.);
  Expects(target_k_eff_tol_ < 1.);
  target_k_eff_lo_ = target_k_eff_ - target_k_eff_tol_;
  target_k_eff_hi_ = target_k_eff_ + target_k_eff_tol_;

  // Get the B10 isotopic abundance to use, if available
  if (node.child("B10_enrichment")) {
    B10_iso_abund_ = node.child("B10_enrichment").text().as_double();
  }
  Expects(B10_iso_abund_ > 0.);
  Expects(B10_iso_abund_ <= 1.);

  // Get the convergence epsilon to use when checking for convergence
  if (node.child("boron_epsilon")) {
    epsilon_ = node.child("boron_epsilon").text().as_double();
  }
  Expects(epsilon_ > 0.);
}

double BoronDriver::solve_ppm(bool first_pass, UncertainDouble& k_eff,
                              UncertainDouble& k_eff_prev)
{
  // The estimation is performed using the Newton-Raphson method
  double del_ppm_del_k;
  double new_ppm;
  if (first_pass) {
    // Make a reasonable update so we can get two points to compute a derivative
    if (k_eff.mean > target_k_eff_) {
      new_ppm = ppm_prev_ + 1000.;
    } else if (k_eff.mean < target_k_eff_) {
      new_ppm = ppm_prev_ - 1000.;
    } else {
      // Then our most recent calculation was on target, so we keep it
      // This should generally never be achieved by a realistic CE MC case due
      // to the noise, but, it may be reached with simpler MG MC cases or
      // deterministic solvers
      new_ppm = ppm_prev_;
    }

  } else {
    del_ppm_del_k = (ppm_ - ppm_prev_) / (k_eff.mean - k_eff_prev.mean);
    new_ppm = (target_k_eff_ - k_eff.mean) * del_ppm_del_k + ppm_;
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

bool BoronDriver::is_converged(UncertainDouble& k_eff) {
  bool ppm_converged;
  bool k_eff_converged;
  double k_eff_95lo;
  double k_eff_95hi;
  std::stringstream msg;
  std::string gt_lt_sign;

  std::ios_base::fmtflags old_flags(std::cout.flags());

  // Check the boron ppm convergence
  ppm_converged = B10_iso_abund_ * (ppm_ - ppm_prev_) < epsilon_;

  // Check to make sure keff is in the expected range to a 95% CI
  k_eff_95lo = k_eff.mean - CI_95 * k_eff.std_dev;
  k_eff_95hi = k_eff.mean + CI_95 * k_eff.std_dev;
  if ((target_k_eff_lo_ <= k_eff_95lo) && (k_eff_95hi <= target_k_eff_hi_)) {
    k_eff_converged = true;
  } else if ((k_eff.mean - target_k_eff_) < target_k_eff_tol_) {
    // If the neutronics k-eff is has higher uncertainty than does the range
    // allows, the above may not work. In that case, lets make sure k-eff is
    // simply in the expected range. If this software were to control the
    // number of histories then this block would have to estimate the new number
    // of histories to simulate to achieve acceptable uncertainties.
    k_eff_converged = true;
  } else {
    k_eff_converged = false;
  }

  msg << "  Boron change: " << std::fixed << std::setprecision(2)
      << B10_iso_abund_ * (ppm_ - ppm_prev_);
  if (ppm_converged) {
    msg << " < ";
  } else {
    msg << " > ";
  }
  msg << epsilon_ << " B10 ppm";
  comm_.message(msg.str());
  msg.clear();
  msg.str(std::string());

  msg << "  k-eff deviation: " << std::fixed << std::setprecision(5)
      << fabs(k_eff.mean - target_k_eff_);
  if (k_eff_converged) {
    msg << " < ";
  } else {
    msg << " > ";
  }
  msg << target_k_eff_tol_;
  comm_.message(msg.str());
  msg.clear();
  msg.str(std::string());

  if (k_eff.std_dev > target_k_eff_tol_) {
    msg << "  k-eff uncertainty (+/-" << k_eff.std_dev << ") > "
        << "target tolerance (" << target_k_eff_tol_ << ")";
    comm_.message(msg.str());
    msg.clear();
    msg.str(std::string());
    comm_.message("  Consider increasing neutronics solver histories");
  }

  std::cout.flags(old_flags);

  return (ppm_converged == true) && (k_eff_converged == true);
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
