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

  // TODO: Implement these
  bool has_coupling_data() const override { return false;}
  int set_heat_source_at(int32_t local_elem, double heat) override {return -1;}
  int n_local_elem() const override {return -1;}
  std::size_t n_global_elem() const override {return -1;};

private:
  std::string setup_file_;
  std::string thread_model_;
  std::string device_number_;
  double time_;
  int tstep_;

  // TODO: Implement these
  std::vector<double> temperature_local() const override {return std::vector<double>{};}
  std::vector<double> density_local() const override {return std::vector<double>{};}
  std::vector<int> fluid_mask_local() const override {return std::vector<int>{};}
  std::vector<Position> centroid_local() const override {return std::vector<Position>{};}
  std::vector<double> volume_local() const override {return std::vector<double>{};}
};

}


#endif // ENRICO_SRC_NEKRS_DRIVER_H

