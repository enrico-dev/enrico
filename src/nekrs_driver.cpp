#include "enrico/nekrs_driver.h"
#include "nekrs.hpp"
#include "gsl.hpp"

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
  Position p{0, 0, 0};
  double mass = 0;
  for (gsl::index i = 0; i < n_gll_; ++i) {
    p.x += nekrs::mesh::x()[local_elem * n_gll_ + i];
    p.y += nekrs::mesh::y()[local_elem * n_gll_ + i];
    p.z += nekrs::mesh::z()[local_elem * n_gll_ + i];
    mass += nekrs::mesh::mm()[local_elem * n_gll_ + i];
  }
  p.x /= mass;
  p.y /= mass;
  p.z /= mass;
  return p;
}

std::vector<Position> NekRSDriver::centroid_local() const {
  std::vector<Position> local_element_centroids(n_local_elem_);
  for (int32_t i = 0; i < n_local_elem_; ++i) {
    local_element_centroids[i] = this->centroid_at(i);
  }
  return local_element_centroids;
}

double NekRSDriver::volume_at(int32_t local_elem) const  {
  Expects(local_elem < n_local_elem_)
  double volume = 0;
  for (int32_t i = 0; i < n_gll_; ++i) {
    volume += nekrs::mesh::mm()[local_elem * n_gll_ + i];
  }
  return volume;
}

std::vector<double> NekRSDriver::volume_local() const {
  std::vector<double> local_elem_volumes(n_local_elem_);
  for (int32_t i = 0; i < n_local_elem_; ++i) {
    local_elem_volumes[i] = this->volume_at(i);
  }
  return local_elem_volumes;
}

}

