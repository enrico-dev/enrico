#include "enrico/nekrs_driver.h"
#include "nekrs.hpp"

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

}
