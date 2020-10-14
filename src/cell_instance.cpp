#include "enrico/cell_instance.h"

#include "enrico/error.h"

#include "openmc/capi.h"

#include <sstream>
#include <stdexcept>

namespace enrico {

CellInstance::CellInstance(Position position)
{
  // Get cell index/instance corresponding to position
  double xyz[3] = {position.x, position.y, position.z};
  err_chk(openmc_find_cell(xyz, &index_, &instance_));

  // Determine what material fills the cell instance
  int type;
  int32_t* indices;
  int32_t n;
  err_chk(openmc_cell_get_fill(index_, &type, &indices, &n));

  // Get index for the appropriate instance
  material_index_ = indices[n > 1 ? instance_ : 0];

  // Get volume of material (if non-void)
  if (material_index_ >= 0) {
    volume_ = this->material()->volume();
  }
}

openmc::Cell* CellInstance::cell() const
{
  return openmc::model::cells.at(index_).get();
}

openmc::Material* CellInstance::material() const
{
  return openmc::model::materials.at(material_index_).get();
}

void CellInstance::set_temperature(double T) const
{
  this->cell()->set_temperature(T, instance_);
}

double CellInstance::get_temperature() const
{
  return this->cell()->temperature(instance_);
}

void CellInstance::set_density(double rho) const
{
  this->material()->set_density(rho, "g/cm3");
}

double CellInstance::get_density() const
{
  return this->material()->density();
}

bool CellInstance::is_fissionable() const
{
  return this->material()->fissionable();
}

bool CellInstance::operator==(const CellInstance& other) const
{
  return index_ == other.index_ && instance_ == other.instance_;
}

CellHandle CellInstance::get_handle() const
{
  // Taken from https://stackoverflow.com/a/17017281
  CellHandle res = 17;
  res = 31 * res + std::hash<int32_t>()(index_);
  res = 31 * res + std::hash<int32_t>()(instance_);
  return res;
}

} // namespace enrico
