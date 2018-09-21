#include "openmc_interface.h"

#include <stdexcept>
#include <sstream>

#include "openmc/capi.h"
#include "error.h"

namespace stream {

CellInstance::CellInstance(Position position)
{
  // Get cell index/instance corresponding to position
  double xyz[3] = {position.x, position.y, position.z};
  err_chk(openmc_find_cell(xyz, &index_, &instance_), openmc_err_msg);

  // Determine what material fills the cell instance
  int type;
  int32_t* indices;
  int32_t n;
  err_chk(openmc_cell_get_fill(index_, &type, &indices, &n), openmc_err_msg);

  // TODO: Right now get_fill returns 0-based indices, but the tally interface
  // expects 1-based. Once tallies move to 0-based, change this.
  material_index_ = indices[instance_] + 1;

  // Get volume of material
  err_chk(openmc_material_get_volume(material_index_, &volume_), openmc_err_msg);
}

openmc::Material* CellInstance::material() const
{
  return openmc::global_materials[material_index_ - 1];
}

void CellInstance::set_temperature(double T) const
{
  openmc_cell_set_temperature(index_, T, &instance_);
}

} // namespace stream
