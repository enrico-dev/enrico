#include "smrt/Assembly_Model.h"

#include "Nemesis/utils/Constants.hh"

namespace stream
{
//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
Assembly_Model::Assembly_Model(const std::vector<PIN_TYPE>& pin_map,
                               const std::vector<double>&   x_edges,
                               const std::vector<double>&   y_edges,
                               double                       height)
    : d_pin_map(pin_map)
    , d_x_edges(x_edges)
    , d_y_edges(y_edges)
    , d_height(height)
{
    d_Nx = d_x_edges.size() - 1;
    d_Ny = d_y_edges.size() - 1;
    Require(d_pin_map.size() == d_Nx * d_Ny);
}

//---------------------------------------------------------------------------//
// Alternate Constructor
//---------------------------------------------------------------------------//
Assembly_Model::Assembly_Model(Teuchos::RCP<Teuchos::ParameterList> params)
{
    using Array_Dbl = Teuchos::Array<double>;
    using Array_Str = Teuchos::Array<std::string>;

    Validate(params->isType<Array_Dbl>("x_edges"),
             "Must specify 'x_edges'.");
    Validate(params->isType<Array_Dbl>("y_edges"),
             "Must specify 'y_edges'.");
    Validate(params->isType<Array_Dbl>("z_edges"),
             "Must specify 'z_edges'.");
    Validate(params->isType<Array_Str>("pin_map"),
             "Must specify 'pin_map'.");

    // Get x/y/z edges
    auto x = params->get<Array_Dbl>("x_edges");
    d_x_edges = x.toVector();
    d_Nx = d_x_edges.size() - 1;
    auto y = params->get<Array_Dbl>("y_edges");
    d_y_edges = y.toVector();
    d_Ny = d_y_edges.size() - 1;
    auto z = params->get<Array_Dbl>("z_edges");
    d_height = z.back();

    // Translate pin map strings to enums
    auto pins = params->get<Array_Str>("pin_map");
    Validate(pins.size() == d_Nx * d_Ny,
             "Expected pin map to be size " << d_Nx * d_Ny
             << " but actual size is " << pins.size());
    d_pin_map.resize(d_Nx*d_Ny);
    for (int pin = 0; pin < d_Nx*d_Ny; ++pin)
    {
        auto val = pins[pin];
        Validate(val == "Fuel" || val == "Guide",
                 "Pin map entries must be 'Fuel' or 'Guide', "
                 " entry " << pin << " is " << val);
        if (val == "Fuel")
            d_pin_map[pin] = FUEL;
        else
            d_pin_map[pin] = GUIDE;
    }

    // Get radii
    d_fuel_radius  = params->get("fuel_radius",  0.418);
    d_clad_radius  = params->get("clad_radius",  0.475);
    d_guide_radius = params->get("guide_radius", 0.612);
}

//---------------------------------------------------------------------------//
// Compute flow area for channel
//---------------------------------------------------------------------------//
double Assembly_Model::flow_area(int i, int j) const
{
    Require(i < d_Nx);
    Require(j < d_Ny);

    double dx = d_x_edges[i+1] - d_x_edges[i];
    double dy = d_y_edges[j+1] - d_y_edges[j];

    auto type = pin_type(i,j);
    double r = 0.0;
    if (type == FUEL)
        r = d_clad_radius;
    else if (type == GUIDE)
        r = d_guide_radius;
    Check(r > 0.0);

    double area = dx * dy - nemesis::constants::pi * r * r;
    Ensure(area > 0);
    return area;
}

//---------------------------------------------------------------------------//
} // end namespace stream

//---------------------------------------------------------------------------//
// end of Assembly_Model.cpp
//---------------------------------------------------------------------------//
