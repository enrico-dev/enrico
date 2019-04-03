#include "enrico/openmc_driver.h"

#include "enrico/const.h"
#include "enrico/error.h"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/tallies/tally.h"
#include "xtensor/xadapt.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xview.hpp"

#include <string>

namespace enrico {

OpenmcDriver::OpenmcDriver(MPI_Comm comm) : Driver(comm)
{
  if (active()) {
    err_chk(openmc_init(0, nullptr, &comm));
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void OpenmcDriver::create_tallies(gsl::span<int32_t> materials)
{
  // Determine maximum tally/filter ID used so far
  int32_t filter_id, tally_id;
  openmc_get_filter_next_id(&filter_id);
  openmc_get_tally_next_id(&tally_id);

  err_chk(openmc_new_filter("material", &index_filter_));
  err_chk(openmc_filter_set_id(index_filter_, filter_id));

  // Set bins for filter
  err_chk(openmc_material_filter_set_bins(index_filter_, materials.size(),
                                          materials.data()));

  // Create tally and assign scores/filters
  int32_t index_tally;
  err_chk(openmc_extend_tallies(1, &index_tally, nullptr));
  err_chk(openmc_tally_set_id(index_tally, tally_id));
  std::vector<std::string> scores {"kappa-fission"};
  tally_ = openmc::model::tallies[index_tally].get();
  tally_->set_scores(scores);
  tally_->set_filters(&index_filter_, 1);
}

xt::xtensor<double, 1> OpenmcDriver::heat_source(double power)
{
  // Determine number of realizatoins for normalizing tallies
  int m = tally_->n_realizations_;

  // Broadcast number of realizations
  // TODO: Change OpenMC so that it's correct on all ranks
  comm_.Bcast(&m, 1, MPI_INT);

  // Determine energy production in each material
  auto mean_value = xt::view(tally_->results_, xt::all(), 0, openmc::RESULT_SUM);
  xt::xtensor<double, 1> heat = JOULE_PER_EV * mean_value / m;

  // Get total heat production [J/source]
  double total_heat = xt::sum(heat)();

  // Normalize heat source in each material and collect in an array
  for (int i = 0; i < cells_.size(); ++i) {
    // Get volume
    double V = cells_.at(i).volume_;

    // Convert heat from [J/source] to [W/cm^3]. Dividing by total_heat gives
    // the fraction of heat deposited in each material. Multiplying by power
    // givens an absolute value in W
    heat(i) *= power / (total_heat * V);
  }

  return heat;
}

void OpenmcDriver::init_step()
{
  err_chk(openmc_simulation_init());
}

void OpenmcDriver::solve_step()
{
  err_chk(openmc_run());
}

void OpenmcDriver::write_step(unsigned int timestep, unsigned int iteration)
{
  std::string filename {"openmc_" + std::to_string(timestep) +
    "_" + std::to_string(iteration) + ".h5"};
  err_chk(openmc_statepoint_write(filename.c_str(), nullptr));
}

void OpenmcDriver::finalize_step()
{
  err_chk(openmc_simulation_finalize());
}

OpenmcDriver::~OpenmcDriver()
{
  if (active()) {
    err_chk(openmc_finalize());
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace enrico
