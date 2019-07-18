//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   shift_driver.cpp
 * \author Steven Hamilton
 * \date   Wed Aug 15 09:25:43 2018
 * \brief  ShiftDriver class definitions.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <map>

#include "enrico/utils.h"
#include "smrt/shift_driver.h"

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Nemesis/comm/global.hh"
#include "Omnibus/driver/Sequence_Shift.hh"
#include "Omnibus/shift_managers/Shift_Tallies.hh"
#include "Shift/mc_tallies/Cell_Union_Tally.hh"

namespace enrico {
//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
ShiftDriver::ShiftDriver(SP_Assembly_Model assembly,
                         std::string shift_input,
                         const std::vector<double>& z_edges)
  : d_assembly(assembly)
  , d_power_tally_name("power")
{
  // Build a Teuchos communicator
  const Teuchos::RCP<const Teuchos::Comm<int>> comm =
    Teuchos::DefaultComm<int>::getComm();

  // Make a temporary Parameter list
  RCP_PL plist = RCP_PL(new Teuchos::ParameterList("Omnibus_plist_root"));

  // Save the input XML path for later output
  plist->set("input_path", shift_input);

  // Load parameters from disk on processor zero and broadcast them
  Teuchos::updateParametersFromXmlFileAndBroadcast(shift_input, plist.ptr(), *comm);

  // Build Problem
  auto problem = std::make_shared<omnibus::Problem>(plist);

  // Build driver
  d_driver = std::make_shared<omnibus::Multiphysics_Driver>(problem);

  // Store geometry
  auto problem_geom = problem->geometry();
  Expects(problem_geom != nullptr);
  d_geometry = std::dynamic_pointer_cast<Geometry>(problem_geom);
  Expects(d_geometry != nullptr);

  // Create or validate power tally
  add_power_tally(plist, z_edges);
}

std::vector<double> ShiftDriver::heat_source(double power) const
{
  // Extract fission rate from Shift tally
  auto sequence = d_driver->sequence();
  auto shift_seq = std::dynamic_pointer_cast<omnibus::Sequence_Shift>(sequence);
  const auto& tallies = shift_seq->tallies();

  std::vector<double> power_by_cell_ID;

  const auto& cell_tallies = tallies.cell_tallies();
  for (const auto& tally : cell_tallies) {
    if (tally->name() == d_power_tally_name) {
      // Tally results are volume-integrated,
      // divide by volume to get volumetric power
      const auto& result = tally->result();
      auto mean = result.mean(0);
      Expects(result.num_multipliers() == 1);

      // Compute global sum for normalization
      double total_power = std::accumulate(mean.begin(), mean.end(), 0.0);

      for (int cellid = 0; cellid < mean.size(); ++cellid) {
        double tally_volume = d_geometry->cell_volume(cellid);
        power_by_cell_ID.push_back(mean[cellid] / tally_volume / total_power);
      }
    }
  }

  double total_power = 0.0;
  for (int cellid = 0; cellid < power_by_cell_ID.size(); ++cellid) {
    total_power += power_by_cell_ID[cellid] * d_geometry->cell_volume(cellid);
  }
  nemesis::global_sum(total_power);

  double norm_factor = power / total_power;
  for (auto& val : power_by_cell_ID)
    val *= norm_factor;

  return power_by_cell_ID;
}

//---------------------------------------------------------------------------//
// Solve
//---------------------------------------------------------------------------//
void ShiftDriver::solve(const std::vector<double>& th_temperature,
                        const std::vector<double>& coolant_density,
                        std::vector<double>& power)
{
  update_temperature(th_temperature);

  // currently does nothing
  update_density(coolant_density);

  // Rebuild problem (loading any new data needed and run transport
  d_driver->rebuild();
  d_driver->run();
}

void ShiftDriver::update_temperature(const std::vector<double>& temperatures)
{
  Expects(temperatures.size() == d_matids.size());
  auto& comps = d_driver->compositions();

  // Average T/H mesh temperatures onto Shift materials
  std::vector<double> material_temperatures(d_num_materials, 0.0);
  for (int elem = 0; elem < temperatures.size(); ++elem) {
    material_temperatures[d_matids[elem]] += d_vfracs[elem] * temperatures[elem];
  }

  nemesis::global_sum(material_temperatures.data(), d_num_materials);

  for (int matid = 0; matid < d_num_materials; ++matid) {
    if (material_temperatures[matid] > 0.0)
      comps[matid]->set_temperature(material_temperatures[matid]);
  }
}
//---------------------------------------------------------------------------//
// Register list of centroids and cell volumes from T/H solver
//---------------------------------------------------------------------------//
void ShiftDriver::set_centroids_and_volumes(
  const std::vector<enrico::Position>& centroids,
  const std::vector<double>& volumes)
{
  assert(centroids.size() == volumes.size());

  d_matids.resize(centroids.size());
  std::vector<int> cells(centroids.size());
  for (int elem = 0; elem < centroids.size(); ++elem) {
    const auto& c = centroids[elem];

    // Find geometric cell and corresponding material
    cells[elem] = d_geometry->find_cell({c.x, c.y, c.z});
    d_matids[elem] = d_geometry->matid(cells[elem]);
  }

  // Compute per-material volume
  d_num_materials = d_driver->compositions().size();
  d_num_shift_cells = d_geometry->num_cells();
  d_power_map.resize(d_num_shift_cells);
  std::vector<double> mat_volumes(d_num_materials, 0.0);
  for (int elem = 0; elem < volumes.size(); ++elem) {
    int matid = d_matids[elem];
    assert(matid < d_num_materials);
    mat_volumes[matid] += volumes[elem];
    d_power_map[cells[elem]].push_back(elem);
  }

  // Perform global reduction of material volumes
  nemesis::global_sum(mat_volumes.data(), d_num_materials);

  // Convert volumes to volume fractions
  d_vfracs.resize(volumes.size());
  for (int elem = 0; elem < volumes.size(); ++elem) {
    int matid = d_matids[elem];
    d_vfracs[elem] = volumes[elem] / mat_volumes[matid];
  }

  int num_cells = d_geometry->num_cells();
  std::vector<double> cell_volumes(num_cells, 0.0);
  for (int elem = 0; elem < centroids.size(); ++elem) {
    const auto& c = centroids[elem];
    int cell = d_geometry->find_cell({c.x, c.y, c.z});
    cell_volumes[cell] += volumes[elem];
  }
}

//---------------------------------------------------------------------------//
// Add power (fission rate) tally to shift problem
//---------------------------------------------------------------------------//
void ShiftDriver::add_power_tally(RCP_PL& pl, const std::vector<double>& z_edges)
{
  Expects(d_assembly != nullptr);
  auto tally_pl = Teuchos::sublist(pl, "TALLY");
  auto cell_pl = Teuchos::sublist(tally_pl, "CELL");

  using Array_Int = Teuchos::Array<int>;
  using Array_Dbl = Teuchos::Array<double>;
  using Array_Str = Teuchos::Array<std::string>;

  if (!cell_pl->isSublist("power")) {
    // If "power" tally doesn't exist, create it
    auto power_pl = Teuchos::sublist(cell_pl, "power");
    power_pl->set("name", d_power_tally_name);
    power_pl->set("normalization", 1.0);
    power_pl->set("description", std::string("power tally"));
    Teuchos::Array<std::string> rxns(1, "fission");
    power_pl->set("reactions", rxns);
    power_pl->set("cycles", std::string("active"));
    Array_Str cells(d_geometry->num_cells());
    for (int cellid = 0; cellid < d_geometry->num_cells(); ++cellid)
      cells[cellid] = std::to_string(cellid);
    Array_Int counts(d_geometry->num_cells(), 1);
    power_pl->set("union_cells", cells);
    power_pl->set("union_lengths", counts);
  } else {
    // If it exists, make sure it aligns with assembly
    auto power_pl = Teuchos::sublist(cell_pl, "power");
    Validate(d_power_tally_name == power_pl->get<std::string>("name"),
             "Incorrect power tally name");
    auto rxns = power_pl->get<Array_Str>("reactions");
    Validate(rxns.size() == 1, "Incorrect number of reactions in power tally");
    Validate(rxns[0] == "fission", "Incorrect reaction in power tally");
    Validate(power_pl->get<std::string>("cycles") == "active",
             "Incorrect cycle designation in power tally.");
    Validate(power_pl->get<std::string>("type") == "grid",
             "Incorrect mesh type in power tally.");

    // Check x edges
    const auto& x_edges = d_assembly->x_edges();
    const auto& x_tally = power_pl->get<Array_Dbl>("x");
    Validate(x_tally.size() == x_edges.size(),
             "Tally specifies incorrect size of x edges");
    for (int i = 0; i < x_edges.size(); ++i)
      Validate(soft_equiv(x_edges[i], x_tally[i]), "Tally specifies incorrect x edge");

    // Check y edges
    const auto& y_edges = d_assembly->y_edges();
    const auto& y_tally = power_pl->get<Array_Dbl>("y");
    Validate(y_tally.size() == y_edges.size(),
             "Tally specifies incorrect size of y edges");
    for (int i = 0; i < y_edges.size(); ++i)
      Validate(soft_equiv(y_edges[i], y_tally[i]), "Tally specifies incorrect y edge");

    // Check z edges
    const auto& z_tally = power_pl->get<Array_Dbl>("z");
    Validate(z_tally.size() == z_edges.size(),
             "Tally specifies incorrect size of z edges");
    for (int i = 0; i < z_edges.size(); ++i)
      Validate(soft_equiv(z_edges[i], z_tally[i]), "Tally specifies incorrect z edge");
  }
}

//---------------------------------------------------------------------------//
} // end namespace enrico

//---------------------------------------------------------------------------//
// end of shift_driver.cpp
//---------------------------------------------------------------------------//
