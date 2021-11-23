#include "enrico/openmc_driver.h"

#include "enrico/const.h"
#include "enrico/error.h"

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/nuclide.h"
#include "openmc/simulation.h"
#include "openmc/summary.h"
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

OpenmcDriver::OpenmcDriver(MPI_Comm comm)
  : NeutronicsDriver(comm)
{
  timer_driver_setup.start();
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
    if (static_cast<openmc::Fill>(type) == openmc::Fill::MATERIAL) {
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

#ifdef _OPENMP
#pragma omp parallel default(none) shared(num_threads)
#pragma omp single
  num_threads = omp_get_num_threads();
#endif

  timer_driver_setup.stop();
}

void OpenmcDriver::create_tallies()
{
  using gsl::index;
  using gsl::narrow_cast;

  // Build vector of material indices on each rank
  // After CoupledDriver::init_mappings, the cells_ array is up-to-date on the root,
  // so we need to send that info to all the other ranks
  std::vector<int32_t> indices;
  std::vector<int32_t> instances;
  if (comm_.is_root()) {
    for (const auto& c : cells_) {
      indices.push_back(c.index_);
      instances.push_back(c.instance_);
    }
  }
  comm_.broadcast(indices);
  comm_.broadcast(instances);
  Ensures(indices.size() == instances.size());

  std::vector<openmc::CellInstance> openmc_instances;
  for (index i = 0; i < indices.size(); ++i) {
    openmc_instances.push_back(
      {narrow_cast<index>(indices[i]), narrow_cast<index>(instances[i])});
  }
  // Create material filter
  auto f = openmc::Filter::create("cellinstance");
  filter_ = dynamic_cast<openmc::CellInstanceFilter*>(f);

  // Set bins for filter
  filter_->set_cell_instances(openmc_instances);

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

  // Determine energy production in each material. Note that xt::view doesn't
  // work with enum
  int i_sum = static_cast<int>(openmc::TallyResult::SUM);
  auto mean_value = xt::view(tally_->results_, xt::all(), 0, i_sum);
  xt::xtensor<double, 1> heat = JOULE_PER_EV * mean_value / m;

  // Get total heat production [J/source]
  double total_heat = xt::sum(heat)();

  // Convert heat from [J/source] to [W/cm^3]
  gsl::index i = 0;
  for (const auto& c : cells_) {
    double V = c.volume_;
    heat.at(i++) *= power / (total_heat * V);
  }
  return heat;
}

std::vector<CellHandle> OpenmcDriver::find(const std::vector<Position>& positions)
{
  std::vector<CellHandle> handles;
  handles.reserve(positions.size());

  for (const auto& r : positions) {
    // Determine cell instance corresponding to global element
    CellInstance c{r};

    // If this cell instance hasn't been saved yet, add it to cells_ and
    // keep track of what index it corresponds to
    auto h = c.get_handle();
    if (cell_index_.find(h) == cell_index_.end()) {
      cell_index_.emplace(h, cells_.size());
      cells_.push_back(c);
    }

    // Set value for cell instance in array
    handles.push_back(h);
  }

  return handles;
}

void OpenmcDriver::set_boron_ppm(double ppm, double H2Odens) const
{
  for (auto& mat : openmc::model::materials) {
    auto nucs = mat->nuclides();
    auto densities = mat->densities();
    double ppmtodens = ppm * H2Odens / (1.0E+6 - ppm);
    double calcdens;

    double N;

    // TODO: Replace 10.8... with a number obtained elsewhere
    N = ppmtodens * openmc::N_AVOGADRO / 10.8118249681472;

    // Is there boron in this material?
    std::vector<std::string> names;
    std::vector<double> new_densities;
    for (int i = 0; i < nucs.size(); i++) {
      int nuc_index = nucs[i];
      auto& nuclide = openmc::data::nuclides[nuc_index];
      double awr = openmc::data::nuclides[nuc_index]->awr_;

      // Add nuclide name to list of names
      names.push_back(nuclide->name_);

      if (nuclide->Z_ == 5) {

        // Calculate density of B10 or B11 corresponding to the given ppm
        // TODO: replace with numbers obtained elsewhere
        if (nuclide->A_ == 10) {
          calcdens = N * 0.1982;
        }
        if (nuclide->A_ == 11) {
          calcdens = N * 0.8018;
        }
        new_densities.push_back(calcdens);
      } else {
        new_densities.push_back(densities[i]);
      }
    }

    mat->set_densities(names, new_densities);
  }
}

void OpenmcDriver::set_density(CellHandle cell, double rho) const
{
  this->cell_instance(cell).material()->set_density(rho, "g/cm3");
}

void OpenmcDriver::set_temperature(CellHandle cell, double T) const
{
  const auto& c = this->cell_instance(cell);
  c.cell()->set_temperature(T, c.instance_);
}

double OpenmcDriver::get_density(CellHandle cell) const
{
  return this->cell_instance(cell).material()->density();
}

double OpenmcDriver::get_temperature(CellHandle cell) const
{
  const auto& c = this->cell_instance(cell);
  return c.cell()->temperature(c.instance_);
}

double OpenmcDriver::get_volume(CellHandle cell) const
{
  return this->cell_instance(cell).volume_;
}

bool OpenmcDriver::is_fissionable(CellHandle cell) const
{
  return this->cell_instance(cell).material()->fissionable();
}

std::string OpenmcDriver::cell_label(CellHandle cell) const
{
  // Get cell instance
  const auto& c = this->cell_instance(cell);

  // Build label
  std::stringstream label;
  label << openmc::model::cells[c.index_]->id_ << " (" << c.instance_ << ")";
  return label.str();
}

gsl::index OpenmcDriver::cell_index(CellHandle cell) const
{
  return cell_index_.at(cell);
}

CellInstance& OpenmcDriver::cell_instance(CellHandle cell)
{
  return cells_.at(cell_index_.at(cell));
}

const CellInstance& OpenmcDriver::cell_instance(CellHandle cell) const
{
  return cells_.at(cell_index_.at(cell));
}

double OpenmcDriver::get_boron_ppm() const
{
  int j = 0;
  int m = 0;
  double ppm = 0;
  double densB10 = 0;
  double densB11 = 0;
  double densB_avg = 0;
  double densH1 = 0;
  double densH2 = 0;
  double densO16 = 0;
  double densO17 = 0;
  double densH2O_avg = 0;

  for (auto& mat : openmc::model::materials) {
    auto nucs = mat->nuclides();
    auto densities = mat->densities();
    auto name = mat->name();
    auto id = mat->id();

    // Is there boron in this material?
    std::vector<std::string> names;

    double B10, B11;
    double H1, H2, O16, O17;

    double test;

    for (int i = 0; i < nucs.size(); i++) {
      int nuc_index = nucs[i];
      auto& nuclide = openmc::data::nuclides[nuc_index];
      double awr = openmc::data::nuclides[nuc_index]->awr_;

      // Add nuclide name to list of names
      names.push_back(nuclide->name_);

      if (nuclide->Z_ == 5) {

        // Calculate density of B10 or B11 corresponding to the given ppm
        if (nuclide->A_ == 10) {

          B10 = densities[i] * (awr)*openmc::MASS_NEUTRON / openmc::N_AVOGADRO;
        }
        if (nuclide->A_ == 11) {
          B11 = densities[i] * (awr)*openmc::MASS_NEUTRON / openmc::N_AVOGADRO;
        }
      }
      if (nuclide->Z_ == 1) {
        if (nuclide->A_ == 1) {
          H1 = densities[i] * (awr)*openmc::MASS_NEUTRON / openmc::N_AVOGADRO;
        }
        if (nuclide->A_ == 2) {
          H2 = densities[i] * (awr)*openmc::MASS_NEUTRON / openmc::N_AVOGADRO;
        }
      }
      if (nuclide->Z_ == 8) {
        if (nuclide->A_ == 16) {
          O16 = densities[i] * (awr)*openmc::MASS_NEUTRON / openmc::N_AVOGADRO;
        }
        if (nuclide->A_ == 17) {

          O17 = densities[i] * (awr)*openmc::MASS_NEUTRON / openmc::N_AVOGADRO;
        }
      }
    }
    // assume elemental boron
    if (O16 > 0 || O17 > 0) {
      double tests = (H1 + H2) / (O16 + O17);
      tests = abs((2.0 - tests) / 2.0) * 100;
      if (tests < 5.0) {
        m = m + 1;
        if (B10 > 0 || B11 > 0) {
          j = j + 1;
        }
        densH1 = densH1 + H1;
        densH2 = densH2 + H2;
        densO16 = densO16 + O16;
        densO17 = densO17 + O17;
        densB10 = densB10 + B10;
        densB11 = densB11 + B11;
      }
      H1 = 0.0, H2 = 0, O16 = 0.0, O17 = 0.0;
      B10 = 0;
      B11 = 0;
      tests = 100.0;
    }
  }

  densB_avg = (densB10 + densB11) / j;
  densH2O_avg = (densH1 + densH2 + densO16 + densO17) / m;
  ppm = densB_avg / (densB_avg + densH2O_avg) * 1.0E+6;
  if (j == 0) {
    ppm = 0;
  };

  return ppm;
}

double OpenmcDriver::get_H2O_dens() const
{
  // TODO Replace iso specific data with data source that isnt hardcoded
  int m = 0;
  double densH1 = 0;
  double densH2 = 0;
  double densO16 = 0;
  double densO17 = 0;
  double densH2O_avg = 0;

  for (auto& mat : openmc::model::materials) {
    auto nucs = mat->nuclides();
    auto densities = mat->densities();
    auto name = mat->name();
    auto id = mat->id();

    double H1, H2, O16, O17;
    double test;

    for (int i = 0; i < nucs.size(); i++) {
      int nuc_index = nucs[i];
      auto& nuclide = openmc::data::nuclides[nuc_index];
      double awr = openmc::data::nuclides[nuc_index]->awr_;

      if (nuclide->Z_ == 1) {
        if (nuclide->A_ == 1) {
          H1 = densities[i] * (awr)*openmc::MASS_NEUTRON / openmc::N_AVOGADRO;
        }
        if (nuclide->A_ == 2) {
          H2 = densities[i] * (awr)*openmc::MASS_NEUTRON / openmc::N_AVOGADRO;
        }
      }
      if (nuclide->Z_ == 8) {
        if (nuclide->A_ == 16) {
          O16 = densities[i] * (awr)*openmc::MASS_NEUTRON / openmc::N_AVOGADRO;
        }
        if (nuclide->A_ == 17) {

          O17 = densities[i] * (awr)*openmc::MASS_NEUTRON / openmc::N_AVOGADRO;
        }
      }
    }
    // assume elemental boron
    if (O16 > 0 || O17 > 0) {
      double tests = (H1 + H2) / (O16 + O17);
      tests = abs((2.0 - tests) / 2.0) * 100;
      if (tests < 5.0) {
        m = m + 1;
        densH1 = densH1 + H1;
        densH2 = densH2 + H2;
        densO16 = densO16 + O16;
        densO17 = densO17 + O17;
      }
      H1 = 0.0, H2 = 0, O16 = 0.0, O17 = 0.0;
      tests = 100.0;
    }
  }

  densH2O_avg = (densH1 + densH2 + densO16 + densO17) / m;

  return densH2O_avg;
}

void OpenmcDriver::init_step()
{
  timer_init_step.start();
  err_chk(openmc_simulation_init());
  timer_init_step.stop();
}

void OpenmcDriver::solve_step()
{
  timer_solve_step.start();
  err_chk(openmc_run());
  err_chk(openmc_reset_timers());
  timer_solve_step.stop();
}

double OpenmcDriver::get_k_effective() const
{
  double keff[2];
  openmc_get_keff(keff);
  return keff[0];
}

void OpenmcDriver::write_step(int timestep, int iteration)
{
  timer_write_step.start();
  std::string suffix{"_t" + std::to_string(timestep) + "_i" + std::to_string(iteration) +
                     ".h5"};
  std::string filename{"openmc" + suffix};
  err_chk(openmc_statepoint_write(filename.c_str(), nullptr));

  std::string prop_file{"properties" + suffix};
  err_chk(openmc_properties_export(prop_file.c_str()));
  timer_write_step.stop();
}

void OpenmcDriver::finalize_step()
{
  timer_finalize_step.start();
  err_chk(openmc_simulation_finalize());
  timer_finalize_step.stop();
}

OpenmcDriver::~OpenmcDriver()
{
  if (active()) {
    err_chk(openmc_finalize());
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace enrico
