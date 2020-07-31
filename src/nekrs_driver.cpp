#include "enrico/nekrs_driver.h"
#include "nekrs.hpp"
#include "gsl.hpp"
#include "iapws/iapws.h"

namespace enrico {
NekRSDriver::NekRSDriver(MPI_Comm comm, pugi::xml_node node) :
  HeatFluidsDriver(comm, node) {

  if (active()) {
    // See vendor/nekRS/src/core/occaDeviceConfig.cpp for valid keys
    setup_file_ = node.child_value("setup_file");
    device_number_ = node.child_value("device_number");
    thread_model_ = node.child_value("thread_model");

    // TODO: Get these values from XML or command line
    int build_only = 0;
    int ci_mode = 0;
    int size_target = 0;
    int debug = 0;

    std::string cache_dir;
    nekrs::setup(comm, build_only, size_target, ci_mode, cache_dir,
                 setup_file_, device_number_,
                 thread_model_);

    // Local and global element counts
    n_local_elem_ = nekrs::mesh::nElements();
    std::size_t n = n_local_elem_;
    MPI_Allreduce(&n, &n_global_elem_, 1, get_mpi_type<std::size_t>(), MPI_SUM, comm_.comm);

    poly_deg_ = nekrs::mesh::polyDeg();
    n_gll_ = (poly_deg_ + 1) + (poly_deg_ + 1) + (poly_deg_ + 1);

    x_ = nekrs::mesh::x();
    y_ = nekrs::mesh::y();
    z_ = nekrs::mesh::z();
    mass_matrix_ = nekrs::mesh::massMatrix();

    // rho energy is field 1 of rho
    rho_energy_ = nekrs::rho()[1 * nekrs::cds::fieldOffset()];

    // Temperature is field 0 of scalarFields
    // TODO: This uses the 1st stage of the time integration.  Is this correct?
    temperature_ = nekrs::cds::scalarFields();
  }
  comm_.Barrier();
}

void NekRSDriver::init_step() {
}

void NekRSDriver::solve_step() {
  const auto start_time = nekrs::startTime();
  const auto final_time = nekrs::finalTime();
  const auto dt = nekrs::dt();

  time_ = start_time;
  tstep_ = 1;

  while ((final_time - time_) / (final_time * dt) > 1e-6){
    nekrs::runStep(time_, dt, tstep_);
    time_ += dt;
    nekrs::udfExecuteStep(time_, tstep_, 0);
    ++tstep_;
  }
}

void NekRSDriver::write_step(int timestep, int iteration) {
  nekrs::copyToNek(timestep, iteration);
  nekrs::nekOutfld();
}

Position NekRSDriver::centroid_at(int32_t local_elem) const {
  Expects(local_elem < n_local_elem_);
  Position c{0., 0., 0.};
  double mass = 0.;
  for (gsl::index i = 0; i < n_gll_; ++i) {
    auto idx = local_elem * n_gll_ + i;
    c.x += x_[idx] * mass_matrix_[idx];
    c.y += y_[idx] * mass_matrix_[idx];
    c.z += z_[idx] * mass_matrix_[idx];
    mass += mass_matrix_[idx];
  }
  c.x /= mass;
  c.y /= mass;
  c.z /= mass;
  return c;
}

std::vector<Position> NekRSDriver::centroid_local() const {
  std::vector<Position> c(n_local_elem_);
  for (int32_t i = 0; i < n_local_elem_; ++i) {
    c[i] = this->centroid_at(i);
  }
  return c;
}

double NekRSDriver::volume_at(int32_t local_elem) const  {
  Expects(local_elem < n_local_elem_)
  double v = 0.;
  for (int32_t i = 0; i < n_gll_; ++i) {
    v += mass_matrix_[local_elem * n_gll_ + i];
  }
  return v;
}

std::vector<double> NekRSDriver::volume_local() const {
  std::vector<double> v(n_local_elem_);
  for (int32_t i = 0; i < n_local_elem_; ++i) {
    v[i] = this->volume_at(i);
  }
  return v;
}

double NekRSDriver::temperature_at(int32_t local_elem) const {
  Expects(local_elem < n_local_elem_);

  double x = 0.;
  double y = 0.;
  for (int32_t i = 0; i < n_gll_; ++i) {
    auto idx = local_elem * n_gll_ + i;
    x += rho_energy_[idx] * mass_matrix_[idx] * temperature_[idx];
    y += rho_energy_[idx] * mass_matrix_[idx];
  }
  return x / y;
}

std::vector<double> NekRSDriver::temperature_local() const {
  std::vector<double> t(n_local_elem_);
  for (int32_t i = 0; i < n_local_elem_; ++i) {
    t[i] = this->temperature_at(i);
  }
  return t;
}


std::vector<double> NekRSDriver::density_local() const
{
  std::vector<double> local_densities(n_local_elem_);

  for (int32_t i = 0; i < n_local_elem_; ++i) {
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


}

