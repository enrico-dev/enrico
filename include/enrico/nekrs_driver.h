#ifndef ENRICO_SRC_NEKRS_DRIVER_H
#define ENRICO_SRC_NEKRS_DRIVER_H

#include "mpi.h"
#include "pugixml.hpp"
#include "enrico/heat_fluids_driver.h"

namespace enrico {
class NekRSDriver : public HeatFluidsDriver {
public:
  NekRSDriver(MPI_Comm comm, pugi::xml_node node);

  // Default destructor is used; NekRS only runs MPI_Finalize
  //~NekRSDriver();

  void init_step() override;
  void solve_step() override;
  void write_step(int timestep, int iteration) override;

  int n_local_elem() const override {return n_local_elem_;}
  std::size_t n_global_elem() const override {return n_global_elem_;};

  Position centroid_at(int32_t local_elem) const;
  std::vector<Position> centroid_local() const override;

  double volume_at(int32_t local_elem) const;
  std::vector<double> volume_local() const override;

  double temperature_at(int32_t local_elem) const;
  std::vector<double> temperature_local() const override;

  std::vector<double> density_local() const override;

  // TODO: Implement these
  bool has_coupling_data() const override { return false;}
  int set_heat_source_at(int32_t local_elem, double heat) override {return -1;}
  int in_fluid_at(int32_t local_elem) const {return 0;};

private:
  std::string setup_file_;
  std::string thread_model_;
  std::string device_number_;
  double time_;
  int tstep_;
  int n_local_elem_;
  std::size_t n_global_elem_;
  int poly_deg_;
  int n_gll_;

  double *x_;
  double *y_;
  double *z_;
  double *mass_matrix_;
  double *temperature_;
  double *rho_energy_;

  // TODO: Implement these
  std::vector<int> fluid_mask_local() const override {return std::vector<int>{};}
};

}


#endif // ENRICO_SRC_NEKRS_DRIVER_H

