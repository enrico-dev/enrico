#include "enrico/openmc_driver.h"

#include "enrico/const.h"
#include "enrico/error.h"

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_material.h"
#include "openmc/tallies/tally.h"
#include "xtensor/xadapt.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xview.hpp"
#include <gsl/gsl>

#include <string>
#include <unordered_map>

namespace enrico {

OpenmcDriver::OpenmcDriver(Comm comm)
  : NeutronicsDriver(comm)
{
  if (active()) {
    err_chk(openmc_init(0, nullptr, &comm));
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // determine number of fissionable cells in model to aid in catching
  // improperly mapped problems
  n_fissionable_cells_ = 0;
  for (gsl::index i = 0; i < openmc::model::cells.size(); ++i) {
    int type;
    int32_t* indices;
    int32_t n;
    err_chk(openmc_cell_get_fill(i, &type, &indices, &n));

    // only check for cells filled with type FILL_MATERIAL (evaluated to '1' enum)
    if (type == openmc::FILL_MATERIAL) {
      for (gsl::index j = 0; j < n; ++j) {
        int material_index = indices[j];

        // skip cells filled with type MATERIAL_VOID (evaluated to '-1' enum)
        if (material_index != -1) {
          const auto& m = openmc::model::materials.at(material_index);

          if (m->fissionable_)
            n_fissionable_cells_++;
        }
      }
    }
  }
}

void OpenmcDriver::create_tallies(std::size_t n)
{
  using gsl::index;
  using gsl::narrow_cast;

  // Build vector of material indices
  std::vector<openmc::CellInstance> instances;
  for (index i = 0; i < n; ++i) {
    const auto& c = cells_[i];
    instances.push_back({narrow_cast<index>(c.index_), narrow_cast<index>(c.instance_)});
  }

  // Create material filter
  auto f = openmc::Filter::create("cellinstance");
  filter_ = dynamic_cast<openmc::CellInstanceFilter*>(f);

  // Set bins for filter
  filter_->set_cell_instances(instances);

  // Create tally and assign scores/filters
  tally_ = openmc::Tally::create();
  tally_->set_scores({"kappa-fission"});
  tally_->add_filter(filter_);
}

xt::xtensor<double, 1> OpenmcDriver::heat_source(double power) const
{
  // Determine number of realizations for normalizing tallies
  int m = tally_->n_realizations_;

  // Broadcast number of realizations
  // TODO: Change OpenMC so that it's correct on all ranks
  comm_.broadcast(m);

  // Determine energy production in each material
  auto mean_value = xt::view(tally_->results_, xt::all(), 0, openmc::RESULT_SUM);
  xt::xtensor<double, 1> heat = JOULE_PER_EV * mean_value / m;

  // Get total heat production [J/source]
  double total_heat = xt::sum(heat)();

  for (gsl::index i = 0; i < heat.size(); ++i) {
    // Get volume
    double V = cells_.at(i).volume_;

    // Convert heat from [J/source] to [W/cm^3]. Dividing by total_heat gives
    // the fraction of heat deposited in each material. Multiplying by power
    // givens an absolute value in W
    heat(i) *= power / (total_heat * V);
  }

  return heat;
}

std::vector<CellHandle> OpenmcDriver::find(const std::vector<Position>& positions)
{
  std::vector<CellHandle> handles;
  handles.reserve(positions.size());

  std::unordered_map<CellInstance, CellHandle> cell_index;

  for (const auto& r : positions) {
    // Determine cell instance corresponding to global element
    CellInstance c{r};

    // If this cell instance hasn't been saved yet, add it to cells_ and
    // keep track of what index it corresponds to
    if (cell_index.find(c) == cell_index.end()) {
      cell_index[c] = cells_.size();
      cells_.push_back(c);
    }
    auto h = cell_index.at(c);

    // Set value for cell instance in array
    handles.push_back(h);
  }
  return handles;
}

void OpenmcDriver::set_density(CellHandle cell, double rho) const
{
  cells_[cell].material()->set_density(rho, "g/cm3");
}

void OpenmcDriver::set_temperature(CellHandle cell, double T) const
{
  const auto& c = cells_[cell];
  c.cell()->set_temperature(T, c.instance_);
}

double OpenmcDriver::get_density(CellHandle cell) const
{
  return cells_[cell].material()->density();
}

double OpenmcDriver::get_temperature(CellHandle cell) const
{
  const auto& c = cells_[cell];
  return c.cell()->temperature(c.instance_);
}

double OpenmcDriver::get_volume(CellHandle cell) const
{
  return cells_[cell].volume_;
}

bool OpenmcDriver::is_fissionable(CellHandle cell) const
{
  return cells_[cell].material()->fissionable();
}

void OpenmcDriver::init_step()
{
  err_chk(openmc_simulation_init());
}

void OpenmcDriver::solve_step()
{
  err_chk(openmc_run());
}

void OpenmcDriver::write_step(int timestep, int iteration)
{
  std::string filename{"openmc_t" + std::to_string(timestep) + "_i" +
                       std::to_string(iteration) + ".h5"};
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
