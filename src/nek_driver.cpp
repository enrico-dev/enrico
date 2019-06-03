#include "enrico/nek_driver.h"

#include "enrico/error.h"
#include "gsl/gsl"
#include "iapws/iapws.h"
#include "nek5000/core/nek_interface.h"

#include <climits>
#include <fstream>
#include <string>
#include <unistd.h>

namespace enrico {

NekDriver::NekDriver(MPI_Comm comm, double pressure_bc, pugi::xml_node node)
  : HeatFluidsDriver(comm, pressure_bc)
{
  lelg_ = nek_get_lelg();
  lelt_ = nek_get_lelt();
  lx1_ = nek_get_lx1();

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

void NekDriver::init_displs()
{
  if (active()) {
    local_counts_.resize(comm_.size);
    local_displs_.resize(comm_.size);

    comm_.Allgather(&nelt_, 1, MPI_INT, local_counts_.data(), 1, MPI_INT);

    local_displs_.at(0) = 0;
    for (int i = 1; i < comm_.size; ++i) {
      local_displs_.at(i) = local_displs_.at(i - 1) + local_counts_.at(i - 1);
    }
  }
}

xt::xtensor<double, 1> NekDriver::temperature() const
{
  // Each Nek proc finds the temperatures of its local elements
  double local_elem_temperatures[nelt_];
  for (int i = 0; i < nelt_; ++i) {
    local_elem_temperatures[i] = this->temperature_at(i + 1);
  }

  xt::xtensor<double, 1> global_elem_temperatures = xt::xtensor<double, 1>();

  // only the rank 0 process allocates the size for the receive buffer
  if (comm_.rank == 0) {
    global_elem_temperatures.resize({gsl::narrow<std::size_t>(nelgt_)});
  }

  // Gather all the local element temperatures onto the root
  comm_.Gatherv(local_elem_temperatures,
                nelt_,
                MPI_DOUBLE,
                global_elem_temperatures.data(),
                local_counts_.data(),
                local_displs_.data(),
                MPI_DOUBLE);

  // only the return value from root should be used, or else a broadcast added here
  return global_elem_temperatures;
}

xt::xtensor<int, 1> NekDriver::fluid_mask() const
{
  int local_fluid_mask[nelt_];
  for (int i = 0; i < nelt_; ++i) {
    local_fluid_mask[i] = this->in_fluid_at(i + 1);
  }

  xt::xtensor<int, 1> global_fluid_mask;
  if (comm_.rank == 0) {
    global_fluid_mask.resize({gsl::narrow<std::size_t>(nelgt_)});
  }

  comm_.Gatherv(local_fluid_mask,
                nelt_,
                MPI_INT,
                global_fluid_mask.data(),
                local_counts_.data(),
                local_displs_.data(),
                MPI_INT);

  return global_fluid_mask;
}

xt::xtensor<double, 1> NekDriver::density() const
{
  double local_densities[nelt_];

  for (int i = 0; i < nelt_; ++i) {
    if (this->in_fluid_at(i + 1) == 1) {
      auto T = this->temperature_at(i + 1);
      // nu1 returns specific volume in [m^3/kg]
      local_densities[i] = 1.0e-3 / iapws::nu1(pressure_bc_, T);
    } else {
      local_densities[i] = 0.0;
    }
  }

  xt::xtensor<double, 1> global_densities;

  if (comm_.rank == 0) {
    global_densities.resize({gsl::narrow<std::size_t>(nelgt_)});
  }

  comm_.Gatherv(local_densities,
                nelt_,
                MPI_DOUBLE,
                global_densities.data(),
                local_counts_.data(),
                local_displs_.data(),
                MPI_DOUBLE);

  return global_densities;
}

void NekDriver::solve_step()
{
  nek_reset_counters();
  C2F_nek_solve();
}

Position NekDriver::centroid_at(int local_elem) const
{
  double x, y, z;
  err_chk(nek_get_local_elem_centroid(local_elem, &x, &y, &z),
          "Could not find centroid of local element " + std::to_string(local_elem));
  return {x, y, z};
}

double NekDriver::volume_at(int local_elem) const
{
  double volume;
  err_chk(nek_get_local_elem_volume(local_elem, &volume),
          "Could not find volume of local element " + std::to_string(local_elem));
  return volume;
}

double NekDriver::temperature_at(int local_elem) const
{
  double temperature;
  err_chk(nek_get_local_elem_temperature(local_elem, &temperature),
          "Could not find temperature of local element " + std::to_string(local_elem));
  return temperature;
}

int NekDriver::in_fluid_at(int local_elem) const
{
  return nek_local_elem_is_in_fluid(local_elem);
}

int NekDriver::set_heat_source_at(int local_elem, double heat)
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
