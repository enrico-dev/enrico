//! \file boron_driver.h
//! Base class for boron criticality search
#ifndef BORON_DRIVER_H
#define BORON_DRIVER_H

#include "enrico/driver.h"
#include "enrico/heat_fluids_driver.h"
#include "enrico/cell_handle.h"

#include <pugixml.hpp>

namespace enrico {
class BoronDriver : public Driver {
public:
  explicit BoronDriver(MPI_Comm comm, pugi::xml_node node);

  ~BoronDriver();

  //! Sets the internal list of cell handles to fluid-bearing cells
  //! \param fluid_cell_handles The handles to the fluid-bearing cells
  void set_fluid_cells(std::vector<CellHandle>& fluid_cell_handles);

  //! Prints the status of boron convergence
  void print_boron();

  //! Estimates the boron concentration in ppm to find criticality condition
  //! \param first_pass If this is the first iteration or not
  //! \param k_eff The latest estimate of k-eff
  //! \param k_eff_prev The previous estimate of k-eff
  //! \return Boron concentration in [ppm]
  double solve_ppm(bool first_pass, UncertainDouble& k_eff,
                   UncertainDouble& k_eff_prev);

  //! Check convergence of the boron concentration
  //! for the current Picard iteration.
  //! \param k_eff The latest estimate of k-eff
  bool is_converged(UncertainDouble& k_eff);

  //! The handles to the fluid cells
  std::vector<CellHandle> fluid_cell_handles_;

  // Isotopic abundance of B10 in Boron
  // The default is the natural abundance from Meija J, Coplen T B, et al,
  // "Isotopic compositions of the elements 2013(IUPAC Technical Report) ",
  // Pure. Appl. Chem. 88 (3), pp.293 - 306(2013).
  // This value is from the best available measurement column and is consistent
  // with the data distriubted with OpenMC
  double B10_iso_abund_ {0.1982};

  // The current and previous values of the boron concentration
  // This is stored as parts-per-million on a number density basis
  double ppm_prev_{0.};
  double ppm_{0.};

  // Denote if this driver is enabled or not
  bool is_enabled_{false};

private:
  double target_k_eff_{1.};
  double target_k_eff_tol_{0.001};
  double epsilon_{1e-3};
  double target_k_eff_lo_{1. - 0.001};
  double target_k_eff_hi_{1. + 0.001};
};
} // namespace enrico

#endif // BORON_DRIVER_H
