#ifndef ENRICO_SRC_NEKRS_DRIVER_H
#define ENRICO_SRC_NEKRS_DRIVER_H

#include "enrico/heat_fluids_driver.h"
#include "mpi.h"
#include "pugixml.hpp"
#include "nrs.hpp"

namespace enrico {
class NekRSDriver : public HeatFluidsDriver {
public:
  NekRSDriver(MPI_Comm comm, pugi::xml_node node);

  ~NekRSDriver();

  void init_step() override;
  void solve_step() override;
  void write_step(int timestep, int iteration) override;

  int n_local_elem() const override { return n_local_elem_; }
  std::size_t n_global_elem() const override { return n_global_elem_; };

  Position centroid_at(int32_t local_elem) const;
  double volume_at(int32_t local_elem) const;
  double temperature_at(int32_t local_elem) const;
  int in_fluid_at(int32_t local_elem) const;

  bool has_coupling_data() const final { return comm_.rank == 0; }

  int set_heat_source_at(int32_t local_elem, double heat) override;

private:
  std::vector<Position> centroid_local() const override;
  std::vector<double> volume_local() const override;
  std::vector<double> temperature_local() const override;
  std::vector<double> density_local() const override;
  std::vector<int> fluid_mask_local() const override;

  void open_lib_udf();
  void close_lib_udf();

  std::string setup_file_;
  std::string thread_model_;
  std::string device_number_;
  double time_;
  int tstep_;
  int n_local_elem_;
  std::size_t n_global_elem_;
  int poly_deg_;
  int n_gll_;

  // TODO: These assume the default values of dfloat, dlong, and hlong in
  //  nekrs/src/libP/include/types.h.  Might want to typedef them
  const double* x_;
  const double* y_;
  const double* z_;
  const double* temperature_;
  const double* rho_cp_;
  const int* element_info_;
  std::vector<double> mass_matrix_;

  nrs_t* nrs_ptr_;

  void* lib_udf_handle_;
  // TODO: Get cache dir from env.  See udfLoadFunction in nekrs/udf/udf.cpp
  const std::string lib_udf_name_ = ".cache/udf/libUDF.so";
  std::vector<double>* localq_;
};

}

#endif // ENRICO_SRC_NEKRS_DRIVER_H
