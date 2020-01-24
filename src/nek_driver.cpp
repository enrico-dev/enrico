#include "enrico/nek_driver.h"

#include "enrico/error.h"
#include "gsl/gsl"
#include "iapws/iapws.h"
#include "nek5000/core/nek_interface.h"
#include "xtensor/xadapt.hpp"

#include <climits>
#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>

namespace enrico {

NekDriver::NekDriver(MPI_Comm comm, double pressure_bc, pugi::xml_node node)
  : HeatFluidsDriver(comm, pressure_bc)
{
  if (active()) {
    casename_ = node.child_value("casename");
    if (comm_.rank == 0) {
      init_session_name();
    }

    MPI_Fint int_comm = MPI_Comm_c2f(comm_.comm);
    C2F_nek_init(static_cast<const int*>(&int_comm));

    nelgt_ = nek_get_nelgt();
    nelt_ = nek_get_nelt();

    init_displs();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void NekDriver::init_session_name()
{
  char path_buffer[PATH_MAX];
  err_chk(getcwd(path_buffer, PATH_MAX) == path_buffer ? 0 : -1,
          "Error writing SESSION.NAME in NekDriver");

  std::ofstream session_name("SESSION.NAME");
  session_name << casename_ << std::endl << path_buffer << std::endl;
  session_name.close();
}

std::vector<double> NekDriver::temperature_local() const
{
  // Each Nek proc finds the temperatures of its local elements
  std::vector<double> local_elem_temperatures(nelt_);
  for (int32_t i = 0; i < nelt_; ++i) {
    local_elem_temperatures[i] = this->temperature_at(i + 1);
  }

  return local_elem_temperatures;
}

std::vector<int> NekDriver::fluid_mask_local() const
{
  std::vector<int> local_fluid_mask(nelt_);
  for (int32_t i = 0; i < nelt_; ++i) {
    local_fluid_mask[i] = this->in_fluid_at(i + 1);
  }
  return local_fluid_mask;
}

std::vector<double> NekDriver::density_local() const
{
  std::vector<double> local_densities(nelt_);

  for (int32_t i = 0; i < nelt_; ++i) {
    if (this->in_fluid_at(i + 1) == 1) {
      auto T = this->temperature_at(i + 1);
      // nu1 returns specific volume in [m^3/kg]
      local_densities[i] = 1.0e-3 / iapws::nu1(pressure_bc_, T);
    } else {
      local_densities[i] = 0.0;
    }
  }

  return local_densities;
}

void NekDriver::solve_step()
{
  nek_reset_counters();
  C2F_nek_solve();
}

Position NekDriver::centroid_at(int32_t local_elem) const
{
  double x, y, z;
  err_chk(nek_get_local_elem_centroid(local_elem, &x, &y, &z),
          "Could not find centroid of local element " + std::to_string(local_elem));
  return {x, y, z};
}

std::vector<Position> NekDriver::centroid_local() const
{
  int n_local = this->n_local_elem();
  std::vector<Position> local_element_centroids(n_local);
  for (int32_t i = 0; i < n_local; ++i) {
    local_element_centroids[i] = this->centroid_at(i + 1);
  }
  return local_element_centroids;
}

double NekDriver::volume_at(int32_t local_elem) const
{
  double volume;
  err_chk(nek_get_local_elem_volume(local_elem, &volume),
          "Could not find volume of local element " + std::to_string(local_elem));
  return volume;
}

std::vector<double> NekDriver::volume_local() const
{
  int n_local = this->n_local_elem();
  std::vector<double> local_elem_volumes(n_local);
  for (int32_t i = 0; i < n_local; ++i) {
    local_elem_volumes[i] = this->volume_at(i + 1);
  }
  return local_elem_volumes;
}

double NekDriver::temperature_at(int32_t local_elem) const
{
  double temperature;
  err_chk(nek_get_local_elem_temperature(local_elem, &temperature),
          "Could not find temperature of local element " + std::to_string(local_elem));
  return temperature;
}

int NekDriver::in_fluid_at(int32_t local_elem) const
{
  return nek_local_elem_is_in_fluid(local_elem);
}

int NekDriver::set_heat_source_at(int32_t local_elem, double heat)
{
  return nek_set_heat_source(local_elem, heat);
}

NekDriver::~NekDriver()
{
  if (active())
    C2F_nek_end();
  MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace enrico
