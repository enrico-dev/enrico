#include "enrico/boron_driver.h"

namespace enrico {

BoronDriver::BoronDriver(MPI_Comm comm)
  : Driver(comm)
{
  MPI_Barrier(MPI_COMM_WORLD);
}

void BoronDriver::set_k_effective(double keff, double keffprev)
{
  k_eff_prev = keffprev;
  k_eff_ = keff;
}

void BoronDriver::set_ppm(double ppm)
{
  ppm_ = ppm;
}

void BoronDriver::set_ppm_prev(double ppmprev)
{
  ppm_prev_ = ppmprev;
}

void BoronDriver::set_H2O_density(double H2Odens)
{
  H2O_dens_ = H2Odens;
}

double BoronDriver::solveppm(int step)
{
  double m;
  m = (ppm_ - ppm_prev_) / (k_eff_ - k_eff_prev);
  if (step == 0) {
    m = -15363.4365923065;
  }
  ppm_ = (1.000 - k_eff_prev) * m + ppm_prev_;

  return ppm_;
}

void BoronDriver::print_boron()
{
  std::ios_base::fmtflags old_flags(std::cout.flags());
  std::stringstream msg;
  msg << "Boron Concentration [ppm]: " << ppm_prev_ << "\n";

  msg << "          Water Density    [g/cm^3]: " << H2O_dens_;
  comm_.message(msg.str());

  std::cout.flags(old_flags);
}
BoronDriver::~BoronDriver()
{

  MPI_Barrier(MPI_COMM_WORLD);
}
} // namespace enrico