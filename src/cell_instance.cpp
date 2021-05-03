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
  // Uses a Canto pairing function:
  // https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
  return (index_ + instance_) * (index_ + instance_ + 1) / 2 + instance_;
}

void CellInstance::invert_handle(CellHandle handle, int32_t& index, int32_t& instance)
{
  // Inverting Cantor pairing function from CellInstance::get_handle()
  // https://en.wikipedia.org/wiki/Pairing_function#Inverting_the_Cantor_pairing_function
  int32_t w = (sqrt(8 * handle + 1) - 1) / 2;
  int32_t t = (w * w + w) / 2;
  instance = handle - t;
  index = w - instance;
}

} // namespace enrico
