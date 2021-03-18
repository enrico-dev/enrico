#include "enrico/shift_driver.h"

#include <gsl/gsl> // for Expects

#include "Omnibus/driver/Sequence_Shift.hh" // for Sequence_Shift
#include "Shift/mc_tallies/Cell_Union_Tally.hh"
#include "Teuchos_DefaultMpiComm.hpp"          // for MpiComm
#include "Teuchos_XMLParameterListHelpers.hpp" // for RCP, ParameterList

#include <unordered_map>

namespace enrico {

ShiftDriver::ShiftDriver(MPI_Comm comm, pugi::xml_node node)
  : NeutronicsDriver{comm}
{
  if (this->active()) {
    // Get Shift filename
    if (!node.child("filename")) {
      throw std::runtime_error{"Must provide Shift filename in enrico.xml"};
    }
    std::string filename = node.child_value("filename");

    // Create a Parameter list; note that it is stored as a member since it is
    // used in the create_tallies method, which takes no arguments
    plist_ = Teuchos::RCP<Teuchos::ParameterList>(
      new Teuchos::ParameterList("Omnibus_plist_root"));

    // Save the input XML path for later output
    plist_->set("input_path", filename);

    // Build a Teuchos communicator
    auto teuchos_comm = Teuchos::MpiComm<int>(comm);

    // Load parameters from disk on processor zero and broadcast them
    Teuchos::updateParametersFromXmlFileAndBroadcast(
      filename, plist_.ptr(), teuchos_comm);

    // Build Problem
    auto problem = std::make_shared<omnibus::Problem>(plist_);

    // Build driver
    driver_ = std::make_shared<omnibus::Multiphysics_Driver>(problem);

    // Store geometry
    auto problem_geom = problem->geometry();
    Expects(problem_geom != nullptr);
    geometry_ = std::dynamic_pointer_cast<geometria::RTK_Core>(problem_geom);
    Expects(geometry_ != nullptr);

    // Initialize number of cells
    num_cells_ = geometry_->num_cells();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////
// NeutronicsDriver interface

std::vector<CellHandle> ShiftDriver::find(const std::vector<Position>& positions)
{
  std::vector<CellHandle> handles;
  handles.reserve(positions.size());

  for (const auto& r : positions) {
    // Find geometric cell and corresponding material
    cell_type cell = geometry_->find_cell({r.x, r.y, r.z});

    // If this cell hasn't been saved yet, add it to cells_ and
    // keep track of what index it corresponds to
    if (cell_index_.find(cell) == cell_index_.end()) {
      cell_index_[cell] = cells_.size();
      cells_.push_back(cell);
    }
    auto h = cell_index_.at(cell);

    // Set value for cell instance in array
    handles.push_back(h);
  }
  return handles;
}

void ShiftDriver::create_tallies()
{
  auto tally_pl = Teuchos::sublist(plist_, "TALLY");
  auto cell_pl = Teuchos::sublist(tally_pl, "CELL");

  using Teuchos::Array;

  auto power_pl = Teuchos::sublist(cell_pl, "power");
  power_pl->set("name", "power");
  power_pl->set("normalization", 1.0);
  power_pl->set("description", std::string("power tally"));
  Array<std::string> rxns(1, "fission");
  power_pl->set("reactions", rxns);
  power_pl->set("cycles", std::string("active"));

  // Create an array of cells -- note that all cells are used here, not just
  // ones that map to a TH element
  Array<std::string> cells(geometry_->num_cells());
  for (int cellid = 0; cellid < geometry_->num_cells(); ++cellid)
    cells[cellid] = std::to_string(cellid);
  Array<int> counts(geometry_->num_cells(), 1);
  power_pl->set("union_cells", cells);
  power_pl->set("union_lengths", counts);
}

void ShiftDriver::set_density(CellHandle handle, double rho) const
{
  Expects(rho > 0);
  auto cell = cells_.at(handle);
  Expects(cell >= 0 && cell < this->n_cells());
  int matid = geometry_->matid(cell);
  driver_->compositions()[matid]->set_density(rho);
}

void ShiftDriver::set_temperature(CellHandle handle, double T) const
{
  Expects(T > 0);
  auto cell = cells_.at(handle);
  Expects(cell >= 0 && cell < this->n_cells());
  int matid = geometry_->matid(cell);
  driver_->compositions()[matid]->set_temperature(T);
}

double ShiftDriver::get_density(CellHandle handle) const
{
  auto cell = cells_.at(handle);
  int matid = geometry_->matid(cell);
  return driver_->compositions()[matid]->density();
}

double ShiftDriver::get_temperature(CellHandle handle) const
{
  auto cell = cells_.at(handle);
  int matid = geometry_->matid(cell);
  return driver_->compositions()[matid]->temperature();
}

double ShiftDriver::get_volume(CellHandle handle) const
{
  auto cell = cells_.at(handle);
  return geometry_->cell_volume(cell);
}

bool ShiftDriver::is_fissionable(CellHandle handle) const
{
  auto cell = cells_.at(handle);
  int matid = geometry_->matid(cell);
  return driver_->compositions()[matid]->is_fissionable();
}

xt::xtensor<double, 1> ShiftDriver::heat_source(double power) const
{
  // Extract fission rate from Shift tally
  auto sequence = driver_->sequence();
  auto shift_seq = std::dynamic_pointer_cast<omnibus::Sequence_Shift>(sequence);
  const auto& tallies = shift_seq->tallies();

  // Initialize an array of zeros for the heat source
  xt::xtensor<double, 1> heat({cells_.size()}, 0.0);

  const auto& cell_tallies = tallies.cell_tallies();
  for (const auto& tally : cell_tallies) {
    if (tally->name() == "power") {
      // Tally results are volume-integrated,
      // divide by volume to get volumetric power
      const auto& result = tally->result();
      auto mean = result.mean(0);
      Expects(result.num_multipliers() == 1);

      // Compute global sum for normalization
      double total_heat = std::accumulate(mean.begin(), mean.end(), 0.0);

      for (cell_type cell = 0; cell < mean.size(); ++cell) {
        // Determine if this cell corresponds to any TH elements
        auto it = cell_index_.find(cell);

        if (it != cell_index_.end()) {
          // Get volume
          double V = geometry_->cell_volume(cell);

          // Get handle for cell
          CellHandle h = it->second;

          // Convert heat to [W/cm^3]. Dividing by total_heat gives the fraction
          // of heat deposited in each material. Multiplying by power gives an
          // absolute value in W.
          heat(h) = power * mean[cell] / (total_heat * V);
        }
      }
    }
    break;
  }

  return heat;
}

std::string ShiftDriver::cell_label(CellHandle handle) const
{
  return std::to_string(cells_.at(handle));
}

////////////////////////////////////////////////////////////////////////////////
// Driver interface

void ShiftDriver::init_step()
{
  timer_init_step.start();
  // Rebuild problem (loading any new data needed and run transport
  driver_->rebuild();
  timer_init_step.stop();
}

void ShiftDriver::solve_step()
{
  timer_solve_step.start();
  driver_->run();
  timer_solve_step.stop();
}

} // end namespace enrico
