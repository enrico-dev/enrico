#include "enrico/shift_driver.h"

#include <gsl/gsl>

namespace enrico {

ShiftDriverNew::ShiftDriverNew(MPI_Comm comm)
  : NeutronicsDriver{comm}
{}

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

}
