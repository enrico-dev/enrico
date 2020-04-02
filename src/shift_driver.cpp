#include "enrico/shift_driver.h"

#include <gsl/gsl> // for Expects

#include "Teuchos_DefaultMpiComm.hpp"          // for MpiComm
#include "Teuchos_XMLParameterListHelpers.hpp" // for RCP, ParameterList

namespace enrico {

ShiftDriverNew::ShiftDriverNew(MPI_Comm comm, pugi::xml_node node)
  : NeutronicsDriver{comm}
{
  // Get Shift filename
  if (!node.child("filename")) {
    throw std::runtime_error{"Must provide Shift filename in enrico.xml"};
  }
  std::string filename = node.child_value("filename");

  // Make a temporary Parameter list
  auto plist = Teuchos::RCP<Teuchos::ParameterList>(
    new Teuchos::ParameterList("Omnibus_plist_root"));

  // Save the input XML path for later output
  plist->set("input_path", filename);

  // Build a Teuchos communicator
  auto teuchos_comm = Teuchos::MpiComm<int>(comm);

  // Load parameters from disk on processor zero and broadcast them
  Teuchos::updateParametersFromXmlFileAndBroadcast(filename, plist.ptr(), teuchos_comm);

  // Build Problem
  auto problem = std::make_shared<omnibus::Problem>(plist);

  // Build driver
  driver_ = std::make_shared<omnibus::Multiphysics_Driver>(problem);

  // Store geometry
  auto problem_geom = problem->geometry();
  Expects(problem_geom != nullptr);
  geometry_ = std::dynamic_pointer_cast<geometria::RTK_Core>(problem_geom);
  Expects(geometry_ != nullptr);

  // Need to initialize:
  //  - matids_
}

////////////////////////////////////////////////////////////////////////////////
// NeutronicsDriver interface

std::vector<CellHandle> ShiftDriverNew::find(const std::vector<Position>& positions)
{
  std::vector<CellHandle> handles;
  handles.reserve(positions.size());
  matids_.reserve(positions.size());

  for (const auto& r : positions) {
    // Find geometric cell and corresponding material
    handles.push_back(geometry_->find_cell({r.x, r.y, r.z}));
    matids_.push_back(geometry_->matid(handles.back()));
  }
  return handles;
}

void ShiftDriverNew::set_density(CellHandle cell, double rho) const
{
  Expects(rho > 0);
  Expects(cell >= 0 && cell < num_shift_cells_);
  int matid = geometry_->matid(cell);
  driver_->compositions()[matid]->set_density(rho);
}

void ShiftDriverNew::set_temperature(CellHandle cell, double T) const
{
  Expects(T > 0);
  Expects(cell >= 0 && cell < num_shift_cells_);
  int matid = geometry_->matid(cell);
  driver_->compositions()[matid]->set_temperature(T);
}

double ShiftDriverNew::get_density(CellHandle cell) const
{
  int matid = geometry_->matid(cell);
  return driver_->compositions()[matid]->density();
}

double ShiftDriverNew::get_temperature(CellHandle cell) const
{
  int matid = geometry_->matid(cell);
  return driver_->compositions()[matid]->temperature();
}

double ShiftDriverNew::get_volume(CellHandle cell) const
{
  return geometry_->cell_volume(cell);
}

bool ShiftDriverNew::is_fissionable(CellHandle cell) const
{
  int matid = geometry_->matid(cell);
  return driver_->compositions()[matid]->is_fissionable();
}

// TODO: Implement
xt::xtensor<double, 1> ShiftDriverNew::heat_source(double power) const {}

// TODO: Implement
void ShiftDriverNew::create_tallies() {}

////////////////////////////////////////////////////////////////////////////////
// Driver interface

// TODO: Implement
void ShiftDriverNew::init_step() {}

// TODO: Implement
void ShiftDriverNew::solve_step() {}

// TODO: Implement
void ShiftDriverNew::write_step(int timestep, int iteration) {}

// TODO: Implement
void ShiftDriverNew::finalize_step() {}

} // end namespace enrico
