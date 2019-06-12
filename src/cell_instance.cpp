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

  // TODO: Right now get_fill returns 0-based indices, but the material API
  // calls expects 1-based. Once those move to 0-based, change this.
  material_index_ = indices[instance_];

  // Get volume of material (if non-void)
  if (material_index_ >= 0) {
    err_chk(openmc_material_get_volume(material_index_, &volume_));
  }
}

openmc::Material* CellInstance::material() const
{
  return openmc::model::materials[material_index_].get();
}

void CellInstance::set_temperature(double T) const
{
  err_chk(openmc_cell_set_temperature(index_, T, &instance_));
}

double CellInstance::get_temperature() const
{
  double T;
  err_chk(openmc_cell_get_temperature(index_, &instance_, &T));
  return T;
}

void CellInstance::set_density(double rho) const
{
  err_chk(openmc_material_set_density(material_index_, rho, "g/cm3"));
}

double CellInstance::get_density() const
{
  double rho;
  err_chk(openmc_material_get_density(material_index));
  return rho;
}

bool CellInstance::operator==(const CellInstance& other) const
{
  return index_ == other.index_ && instance_ == other.instance_;
}

} // namespace enrico
