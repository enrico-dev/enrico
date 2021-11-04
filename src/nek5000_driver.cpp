#include "enrico/nek5000_driver.h"

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

Nek5000Driver::Nek5000Driver(MPI_Comm comm, pugi::xml_node node)
  : HeatFluidsDriver(comm, node)
{
  timer_driver_setup.start();
  if (active()) {
    casename_ = node.child_value("casename");
    if (node.child("output_heat_source")) {
      output_heat_source_ = node.child("output_heat_source").text().as_bool();
    }

    if (comm_.rank == 0) {
      init_session_name();
    }

    MPI_Fint int_comm = MPI_Comm_c2f(comm_.comm);
    C2F_nek_init(static_cast<const int*>(&int_comm));

    nelgt_ = nek_get_nelgt();
    nelt_ = nek_get_nelt();
    ldimt_ = nek_get_ldimt();
    npscal_ = nek_get_npscal();

    // Check that we have enough storage space for localq.  We store
    // localq as the last passive scalar (hence ldimt >= 2) but we must
    // ensure that it is not used in the solver (hence npascl < ldimt - 1)
    if (ldimt_ < 2 || npscal_ >= ldimt_ - 1) {
      std::stringstream msg;
      msg << "User specified ldimt=" << ldimt_ << " and npscal=" << npscal_
          << ".  For coupling, ENRICO requires ldimt >= 2 and npscal < ldimt - 1";
      std::runtime_error(msg.str());
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  timer_driver_setup.stop();
}

void Nek5000Driver::init_session_name()
{
  char path_buffer[PATH_MAX];
  err_chk(getcwd(path_buffer, PATH_MAX) == path_buffer ? 0 : -1,
          "Error writing SESSION.NAME in NekDriver");

  std::ofstream session_name("SESSION.NAME");
  session_name << casename_ << std::endl << path_buffer << std::endl;
  session_name.close();
}

std::vector<double> Nek5000Driver::temperature() const
{
  // Each Nek proc finds the temperatures of its local elements
  std::vector<double> local_elem_temperatures(nelt_);
  for (int32_t i = 0; i < nelt_; ++i) {
    local_elem_temperatures[i] = this->temperature_at(i);
  }

  return local_elem_temperatures;
}

std::vector<int> Nek5000Driver::fluid_mask() const
{
  std::vector<int> local_fluid_mask(nelt_);
  for (int32_t i = 0; i < nelt_; ++i) {
    local_fluid_mask[i] = this->in_fluid_at(i);
  }
  return local_fluid_mask;
}

std::vector<double> Nek5000Driver::density() const
{
  std::vector<double> local_densities(nelt_);

  for (int32_t i = 0; i < nelt_; ++i) {
    if (this->in_fluid_at(i) == 1) {
      auto T = this->temperature_at(i);
      // nu1 returns specific volume in [m^3/kg]
      local_densities[i] = 1.0e-3 / iapws::nu1(pressure_bc_, T);
    } else {
      local_densities[i] = 0.0;
    }
  }

  return local_densities;
}

void Nek5000Driver::solve_step()
{
  timer_solve_step.start();
  nek_reset_counters();
  C2F_nek_solve();
  timer_solve_step.stop();
}

Position Nek5000Driver::centroid_at(int32_t local_elem) const
{
  double x, y, z;
  err_chk(nek_get_local_elem_centroid(local_elem + 1, &x, &y, &z),
          "Could not find centroid of local element " + std::to_string(local_elem));
  return {x, y, z};
}

std::vector<Position> Nek5000Driver::centroid() const
{
  int n_local = this->n_local_elem();
  std::vector<Position> local_element_centroids(n_local);
  for (int32_t i = 0; i < n_local; ++i) {
    local_element_centroids[i] = this->centroid_at(i);
  }
  return local_element_centroids;
}

double Nek5000Driver::volume_at(int32_t local_elem) const
{
  double volume;
  err_chk(nek_get_local_elem_volume(local_elem + 1, &volume),
          "Could not find volume of local element " + std::to_string(local_elem));
  return volume;
}

std::vector<double> Nek5000Driver::volume() const
{
  int n_local = this->n_local_elem();
  std::vector<double> local_elem_volumes(n_local);
  for (int32_t i = 0; i < n_local; ++i) {
    local_elem_volumes[i] = this->volume_at(i);
  }
  return local_elem_volumes;
}

double Nek5000Driver::temperature_at(int32_t local_elem) const
{
  double temperature;
  err_chk(nek_get_local_elem_temperature(local_elem + 1, &temperature),
          "Could not find temperature of local element " + std::to_string(local_elem));
  return temperature;
}

int Nek5000Driver::in_fluid_at(int32_t local_elem) const
{
  return nek_local_elem_is_in_fluid(local_elem + 1);
}

int Nek5000Driver::set_heat_source_at(int32_t local_elem, double heat)
{
  Expects(local_elem >= 0 && local_elem < nelt_);
  return nek_set_heat_source(local_elem + 1, heat);
}

void Nek5000Driver::write_step(int timestep, int iteration)
{
  nek_write_step(int(output_heat_source_));
}

Nek5000Driver::~Nek5000Driver()
{
  if (active())
    C2F_nek_end();
  MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace enrico
