#include "stream/openmc_driver.h"

#include "stream/const.h"
#include "stream/error.h"

#include "openmc/capi.h"
#include "xtensor/xadapt.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xview.hpp"

namespace stream {

OpenmcDriver::OpenmcDriver(MPI_Comm comm) : TransportDriver(comm)
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

  err_chk(openmc_extend_filters(1, &index_filter_, nullptr));
  err_chk(openmc_filter_set_type(index_filter_, "material"));
  err_chk(openmc_filter_set_id(index_filter_, filter_id));

  // Set bins for filter
  err_chk(openmc_material_filter_set_bins(index_filter_, materials.size(),
                                          materials.data()));

  // Create tally and assign scores/filters
  err_chk(openmc_extend_tallies(1, &index_tally_, nullptr));
  err_chk(openmc_tally_allocate(index_tally_, "generic"));
  err_chk(openmc_tally_set_id(index_tally_, tally_id));
  char score_array[][20]{"kappa-fission"};
  const char* scores[]{score_array[0]}; // OpenMC expects a const char**, ugh
  err_chk(openmc_tally_set_scores(index_tally_, 1, scores));
  err_chk(openmc_tally_set_filters(index_tally_, 1, &index_filter_));
}

xt::xtensor<double, 3> OpenmcDriver::tally_results()
{
  // Get material bins
  int32_t* mats;
  int32_t n_mats;
  err_chk(openmc_material_filter_get_bins(index_filter_, &mats, &n_mats));

  // Get tally results and number of realizations
  double* results;
  int shape_int[3];
  err_chk(openmc_tally_results(index_tally_, &results, shape_int));
  int32_t m;
  err_chk(openmc_tally_get_n_realizations(index_tally_, &m));

  // Determine shape and size
  // TODO: Change the order of shape in OpenMC itself so we don't have to reverse it here
  std::vector<int> shape {shape_int[2], shape_int[1], shape_int[0]};
  int size {shape_int[0] * shape_int[1] * shape_int[2]};

  // Adapt array into xtensor with no ownership
  return xt::adapt(results, size, xt::no_ownership(), shape);
}

xt::xtensor<double, 1> OpenmcDriver::heat_source(double power)
{
  // Get tally results
  auto results {tally_results()};

  // Determine number of realizatoins for normalizing tallies
  int32_t m;
  err_chk(openmc_tally_get_n_realizations(index_tally_, &m));

  // Determine energy production in each material
  auto mean_value = xt::view(results, xt::all(), 0, 1);
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
  // TODO: OpenMC should properly reset tallies/realizations when initializing a
  // simulation
  err_chk(openmc_reset());
}

void OpenmcDriver::solve_step() { err_chk(openmc_run()); }

void OpenmcDriver::finalize_step() { err_chk(openmc_simulation_finalize()); }

OpenmcDriver::~OpenmcDriver()
{
  if (active()) {
    err_chk(openmc_finalize());
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace stream
