//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Shift_Solver.cpp
 * \author Steven Hamilton
 * \date   Wed Aug 15 09:25:43 2018
 * \brief  Shift_Solver class definitions.
 * \note   Copyright (c) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "smrt/Shift_Solver.h"

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Shift/mc_tallies/Cartesian_Mesh_Tally.hh"
#include "Omnibus/driver/Sequence_Shift.hh"
#include "Omnibus/shift_managers/Shift_Tallies.hh"

namespace stream
{
//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
Shift_Solver::Shift_Solver(SP_Assembly_Model          assembly,
                           std::string                shift_input,
                           const std::vector<double>& z_edges)
    : d_assembly(assembly)
    , d_power_tally_name("power")
{
    // Build a Teuchos communicator
    const Teuchos::RCP<const Teuchos::Comm<int> > comm
        = Teuchos::DefaultComm<int>::getComm();

    // Make a temporary Parameter list
    RCP_PL plist = RCP_PL(new Teuchos::ParameterList("Omnibus_plist_root"));

    // Save the input XML path for later output
    plist->set("input_path", shift_input);

    // Load parameters from disk on processor zero and broadcast them
    Teuchos::updateParametersFromXmlFileAndBroadcast(
            shift_input, plist.ptr(), *comm);

    add_power_tally(plist, z_edges);

    // Build Problem
    auto problem = std::make_shared<omnibus::Problem>(plist);

    // Build driver
    d_driver = std::make_shared<omnibus::Multiphysics_Driver>(problem);

    // Store geometry
    auto problem_geom = problem->geometry();
    Check(problem_geom);
    d_geometry = std::dynamic_pointer_cast<Geometry>(problem_geom);
    Check(d_geometry);
}

//---------------------------------------------------------------------------//
// Solve
//---------------------------------------------------------------------------//
void Shift_Solver::solve(
        const std::vector<double>& fuel_temperature,
        const std::vector<double>& coolant_density,
              std::vector<double>& power)
{
    auto& comps = d_driver->compositions();

    // Loop through geometry and update fuel temperature and moderator density
    const auto& core_array = d_geometry->array();
    for (int k = 0; k < core_array.size(def::K); ++k)
    {
        for (int j = 0; j < core_array.size(def::J); ++j)
        {
            for (int i = 0; i < core_array.size(def::I); ++i)
            {
                const auto& assbly_array = core_array.object(i,j,k);
                for (int kk = 0; kk < assbly_array.size(def::K); ++kk)
                {
                    for (int jj = 0; jj < assbly_array.size(def::J); ++jj)
                    {
                        for (int ii = 0; ii < assbly_array.size(def::I); ++ii)
                        {
                            auto type = d_assembly->pin_type(ii,jj);

                            // NOTE: By convention, RTK assemblies are built as
                            // 2-D objects (1 axial level) with axial variation
                            // handled at the core level.  Therefore, the RTK
                            // object using the assembly axial index while the
                            // SMRT model uses the core-level index.
                            const auto& pin = assbly_array.object(ii, jj, kk);
                            int region_id = d_assembly->index(ii, jj, k);

                            if (type == Assembly_Model::FUEL)
                            {
                                // Change fuel temperature
                                int fuel_mat = pin.matid(0);
                                Check(fuel_mat < comps.size());
                                comps[fuel_mat]->set_temperature(
                                    fuel_temperature[region_id]);

                                // Change coolant density
                                int mod_mat  = pin.matid(pin.num_regions()-1);
                                Check(mod_mat < comps.size());
                                comps[mod_mat]->set_density(
                                    coolant_density[region_id]);
                            }
                            else
                            {
                                // Change coolant density inside guide tube
                                int inner_mod_mat  = pin.matid(0);
                                Check(inner_mod_mat < comps.size());
                                comps[inner_mod_mat]->set_density(
                                    coolant_density[region_id]);

                                // Change coolant density outside guide tube
                                int outer_mod_mat  = pin.matid(
                                    pin.num_regions()-1);
                                Check(outer_mod_mat < comps.size());
                                comps[outer_mod_mat]->set_density(
                                    coolant_density[region_id]);
                            }
                        }
                    }
                }
            }
        }
    }

    // Rebuild problem (loading any new data needed and run transport
    d_driver->rebuild();
    d_driver->run();

    // Extract fission rate from Shift tally
    auto sequence = d_driver->sequence();
    auto shift_seq =
        std::dynamic_pointer_cast<omnibus::Sequence_Shift>(sequence);
    const auto& tallies = shift_seq->tallies();

    const auto& mesh_tallies = tallies.mesh_tallies();
    for (const auto& tally : mesh_tallies)
    {
        if (tally->name() == d_power_tally_name)
        {
            // Tally results are volume-integrated,
            // no need to weight with volume
            const auto& result = tally->result();
            auto mean = result.mean(0);
            Check(result.num_multipliers() == 1);
            Check(mean.size() == power.size());
            std::copy(mean.begin(), mean.end(), power.begin());
        }
    }
}

//---------------------------------------------------------------------------//
// Add power (fission rate) tally to shift problem
//---------------------------------------------------------------------------//
void Shift_Solver::add_power_tally(
        RCP_PL&                    pl,
        const std::vector<double>& z_edges)
{
    Require(d_assembly);
    auto tally_pl = Teuchos::sublist(pl, "TALLY");
    auto mesh_pl = Teuchos::sublist(tally_pl, "MESH");

    using Array_Dbl = Teuchos::Array<double>;
    using Array_Str = Teuchos::Array<std::string>;

    if (!mesh_pl->isSublist("power"))
    {
        // If "power" tally doesn't exist, create it
        auto power_pl = Teuchos::sublist(mesh_pl, "power");
        power_pl->set("name", d_power_tally_name);
        power_pl->set("normalization", 1.0);
        power_pl->set("description", std::string("power tally"));
        Teuchos::Array<std::string> rxns(1, "fission");
        power_pl->set("reactions", rxns);
        power_pl->set("cycles", std::string("active"));
        power_pl->set("type", std::string("grid"));
        Array_Dbl x(d_assembly->x_edges().begin(),
                    d_assembly->x_edges().end());
        power_pl->set("x", x);
        Array_Dbl y(d_assembly->y_edges().begin(),
                    d_assembly->y_edges().end());
        power_pl->set("y", y);
        Array_Dbl z(z_edges.begin(), z_edges.end());
        power_pl->set("z", z);
    }
    else
    {
        // If it exists, make sure it aligns with assembly
        auto power_pl = Teuchos::sublist(mesh_pl, "power");
        Validate(d_power_tally_name == power_pl->get<std::string>("name"),
                 "Incorrect power tally name");
        auto rxns = power_pl->get<Array_Str>("reactions");
        Validate(rxns.size() == 1,
                 "Incorrect number of reactions in power tally");
        Validate(rxns[0] == "fission",
                 "Incorrect reaction in power tally");
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
            Validate(nemesis::soft_equiv(x_edges[i], x_tally[i]),
                     "Tally specifies incorrect x edge");

        // Check y edges
        const auto& y_edges = d_assembly->y_edges();
        const auto& y_tally = power_pl->get<Array_Dbl>("y");
        Validate(y_tally.size() == y_edges.size(),
                 "Tally specifies incorrect size of y edges");
        for (int i = 0; i < y_edges.size(); ++i)
            Validate(nemesis::soft_equiv(y_edges[i], y_tally[i]),
                     "Tally specifies incorrect y edge");

        // Check z edges
        const auto& z_tally = power_pl->get<Array_Dbl>("z");
        Validate(z_tally.size() == z_edges.size(),
                 "Tally specifies incorrect size of z edges");
        for (int i = 0; i < z_edges.size(); ++i)
            Validate(nemesis::soft_equiv(z_edges[i], z_tally[i]),
                     "Tally specifies incorrect z edge");
    }
}

//---------------------------------------------------------------------------//
} // end namespace stream

//---------------------------------------------------------------------------//
// end of Shift_Solver.cpp
//---------------------------------------------------------------------------//
