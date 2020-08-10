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
  double volume_at(int32_t local_elem) const;
  double temperature_at(int32_t local_elem) const;
  int in_fluid_at(int32_t local_elem) const;

  bool has_coupling_data() const final { return comm_.rank == 0; }

  // TODO: Implement this
  int set_heat_source_at(int32_t local_elem, double heat) override {return -1;}

private:
  std::vector<Position> centroid_local() const override;
  std::vector<double> volume_local() const override;
  std::vector<double> temperature_local() const override;
  std::vector<double> density_local() const override;
  std::vector<int> fluid_mask_local() const override;

  std::string setup_file_;
  std::string thread_model_;
  std::string device_number_;
  double time_;
  int tstep_;
  int n_local_elem_;
  std::size_t n_global_elem_;
  int poly_deg_;
  int n_gll_;

  const double* x_;
  const double* y_;
  const double* z_;
  const double* mass_matrix_;
  const double* temperature_;
  const double* rho_energy_;
  const int* element_info_;
};

}


#endif // ENRICO_SRC_NEKRS_DRIVER_H

