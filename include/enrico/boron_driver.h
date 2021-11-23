//! \file boron_driver.h
//! Base class for boron criticality search
#ifndef BORON_DRIVER_H
#define BORON_DRIVER_H

#include "enrico/driver.h"

namespace enrico {
class BoronDriver : public Driver {
public:
  explicit BoronDriver(MPI_Comm comm);

  ~BoronDriver();

  //! Gets the water density of the borated water
  //! \return Water density in [g/cm^3]
  void set_H2O_density(double H2Odens);

  //! Sets the current and previous keff in the boron driver class
  //! \param keff is the current k-effective after the Openmc run
  //! \param keffprev is the previous k-effective
  void set_k_effective(double keff, double keffprev);

  //! Sets the current boron concentration in ppm in the boron driver Openmc class
  //! \param ppm is the current boron concentration in [ppm] after the Openmc run
  void set_ppm(double ppm);

  //! Sets the previous boron concentration in ppm in the boron driver Openmc class
  //! \param ppm is the previous boron concentration in [ppm] in the previous Openmc run
  void set_ppm_prev(double ppmprev);

  //! Prints the current boron concentration in [ppm] and the current water density in
  //! [g/cm^3]
  void print_boron();

  //! Estimates the boron concentration in ppm to find criticality condition
  //! \param step is used to determine if an initial slope is needed
  //! \return Boron concentration in [ppm]
  double solveppm(int step);

private:
  double k_eff_;
  double k_eff_prev;
  double ppm_;
  double ppm_prev_;
  double H2O_dens_;
};
} // namespace enrico

#endif // BORON_DRIVER_H