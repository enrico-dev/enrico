#include "enrico/nekrs_driver.h"
#include "enrico/error.h"
#include "gsl.hpp"
#include "iapws/iapws.h"
#include "libP/include/mesh3D.h"
#include "nekrs.hpp"

#include <dlfcn.h>
#include <algorithm>

namespace enrico {
NekRSDriver::NekRSDriver(MPI_Comm comm, pugi::xml_node node)
  : HeatFluidsDriver(comm, node)
{

  if (active()) {
    // See vendor/nekRS/src/core/occaDeviceConfig.cpp for valid keys
    setup_file_ = node.child_value("casename");
    nekrs::setup(comm, 0, 0, 0, "", setup_file_, "", "");

    open_lib_udf();

    // Check that we're running a CHT simulation.
    err_chk(nekrs::cht(), "NekRS simulation is not setup for conjugate-heat transfer (CHT).  "
            "ENRICO must be run with a CHT simulation.");

    // Local and global element counts
    n_local_elem_ = nekrs::mesh::nElements();
    std::size_t n = n_local_elem_;
    MPI_Allreduce(
      &n, &n_global_elem_, 1, get_mpi_type<std::size_t>(), MPI_SUM, comm_.comm);

    poly_deg_ = nekrs::mesh::polyDeg();
    n_gll_ = (poly_deg_ + 1) * (poly_deg_ + 1) * (poly_deg_ + 1);

    x_ = nekrs::mesh::x();
    y_ = nekrs::mesh::y();
    z_ = nekrs::mesh::z();
    element_info_ = nekrs::mesh::elementInfo();

    // rho energy is field 1 (0-based) of rho
    rho_cp_ = nekrs::rhoCp();

    temperature_ = nekrs::scalarFields();

    // Construct lumped mass matrix from vgeo
    // See cdsSetup in vendor/nekRS/src/core/insSetup.cpp
    mass_matrix_.resize(n_local_elem_ * n_gll_);
    for (gsl::index e = 0; e < n_local_elem_; ++e) {
      for (gsl::index n = 0; n < n_gll_; ++n) {
        mass_matrix_[e * n_gll_ + n] =
          nekrs::mesh::vgeo()[e * n_gll_ * nekrs::mesh::nVgeo() + JWID * n_gll_ + n];
      }
    }

    // for (gsl::index i = 0; i < mass_matrix_.size(); ++i) {
    //  std::cout << "Rank, elem, mass: " << comm_.rank << ", " << i << ", " <<
    //  mass_matrix_[i] << std::endl; if (mass_matrix_[i] < 1e-16) {
    //    std::runtime_error("Very small mass at rank, elem: " +
    //                       std::to_string(comm_.rank) + ", " + std::to_string(i));
    //  }
    //}
    comm_.message("size of local_q: " + std::to_string(localq_->size()));
    comm_.message("N local elem: " + std::to_string(n_local_elem_));

    init_displs();
  }
  comm_.Barrier();
}

void NekRSDriver::init_step() {
  auto min = std::min_element(localq_->cbegin(), localq_->cend());
  auto max = std::max_element(localq_->cbegin(), localq_->cend());
  std::cout << "[ENRICO] : Min, max localq at " << __FILE__ << ":" << __LINE__ << ": " 
    << *min << ", " << *max << std::endl;
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
  // LOOK HERE!!!
  nekrs::copyToNek(time_, tstep_);
}

void NekRSDriver::write_step(int timestep, int iteration) {
  nekrs::copyToNek(timestep, iteration);
  nekrs::nekOutfld();
  return; 
}

Position NekRSDriver::centroid_at(int32_t local_elem) const {
  Expects(local_elem < n_local_elem());
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
  std::vector<Position> c(n_local_elem());
  for (int32_t i = 0; i < n_local_elem(); ++i) {
    c[i] = this->centroid_at(i);
  }
  return c;
}

double NekRSDriver::volume_at(int32_t local_elem) const
{
  Expects(local_elem < n_local_elem());
  double v = 0.;
  for (int32_t i = 0; i < n_gll_; ++i) {
    v += mass_matrix_[local_elem * n_gll_ + i];
  }
  return v;
}

std::vector<double> NekRSDriver::volume_local() const {
  std::vector<double> v(n_local_elem());
  for (int32_t i = 0; i < n_local_elem(); ++i) {
    v[i] = this->volume_at(i);
  }
  return v;
}

double NekRSDriver::temperature_at(int32_t local_elem) const {
  Expects(local_elem < n_local_elem());

  double sum0 = 0.;
  double sum1 = 0.;
  for (int32_t i = 0; i < n_gll_; ++i) {
    auto idx = local_elem * n_gll_ + i;
    sum0 += rho_cp_[idx] * temperature_[idx];
    sum1 += rho_cp_[idx];
  }
  return sum0 / sum1;
}

std::vector<double> NekRSDriver::temperature_local() const 
{
  std::vector<double> t(n_local_elem());
  for (int32_t i = 0; i < n_local_elem(); ++i) {
    t[i] = this->temperature_at(i);
  }

  auto min = std::min_element(t.cbegin(), t.cend());
  auto max = std::max_element(t.cbegin(), t.cend());
  std::cout << "[ENRICO] : Min, max t at " << __FILE__ << ":" << __LINE__ << ": " 
    << *min << ", " << *max << std::endl;

  return t;
}

std::vector<double> NekRSDriver::density_local() const
{
  nekrs::copyToNek(time_, tstep_);
  std::vector<double> local_densities(n_local_elem());

  for (int32_t i = 0; i < n_local_elem(); ++i) {
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

int NekRSDriver::in_fluid_at(int32_t local_elem) const {
  // In NekRS, element_info_[i] == 1 if i is a *solid* element
  Expects(local_elem < n_local_elem());
  return element_info_[local_elem] == 1 ? 0 : 1;
}

std::vector<int> NekRSDriver::fluid_mask_local() const
{
  std::vector<int> mask(n_local_elem());
  for (int32_t i = 0; i < n_local_elem(); ++i) {
    mask[i] = in_fluid_at(i);
  }
  return mask;
}

int NekRSDriver::set_heat_source_at(int32_t local_elem, double heat)
{
  // TODO: When PR #111 is approved, we won't need this one-off adjustment
  auto e = local_elem - 1;
  for (int i = 0; i < n_gll_; ++i) {
    localq_->at(e * n_gll_ + i) = heat;
  }
  return 0;
}

void NekRSDriver::open_lib_udf()
{
  lib_udf_handle_ = dlopen(lib_udf_name_.c_str(), RTLD_LAZY);
  if (lib_udf_handle_ == NULL) {
    throw std::runtime_error(std::string{"dlopen error for localq: "} + dlerror());
  }
  void* localq_void = dlsym(lib_udf_handle_, "localq");
  if (dlerror() != NULL) {
    throw std::runtime_error("dlsym error for localq");
  }
  localq_ = reinterpret_cast<std::vector<double>*>(localq_void);
}

void NekRSDriver::close_lib_udf()
{
  dlclose(lib_udf_handle_);
}

NekRSDriver::~NekRSDriver()
{
  close_lib_udf();
}

}
