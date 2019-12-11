#include "enrico/openmc_heat_driver.h"

#include "enrico/error.h"
#include "enrico/message_passing.h"
#include "gsl/gsl"
#include "openmc/constants.h"
#include "openmc/tallies/filter_cell_instance.h"
#include "xtensor/xstrided_view.hpp"

#include <cmath>
#include <unordered_map>

namespace enrico {

OpenmcHeatDriver::OpenmcHeatDriver(MPI_Comm comm, pugi::xml_node node)
  : CoupledDriver(comm, node)
{
  // Initialize OpenMC and surrogate heat drivers
  openmc_driver_ = std::make_unique<OpenmcDriver>(comm);
  pugi::xml_node surr_node = node.child("heat_surrogate");
  double pressure_bc = node.child("pressure_bc").text().as_double();
  heat_driver_ = std::make_unique<SurrogateHeatDriver>(comm, pressure_bc, surr_node);

  // Create mappings for fuel pins and setup tallies for OpenMC
  init_mappings();
  init_tallies();

  init_temperatures();
  init_densities();
  init_heat_source();
}

bool OpenmcHeatDriver::has_global_coupling_data() const
{
  return openmc_driver_->active() && heat_driver_->active();
}

NeutronicsDriver& OpenmcHeatDriver::get_neutronics_driver() const
{
  return *openmc_driver_;
}

HeatFluidsDriver& OpenmcHeatDriver::get_heat_driver() const
{
  return *heat_driver_;
}

void OpenmcHeatDriver::init_mappings()
{
  using openmc::PI;

  const auto& r_fuel = heat_driver_->r_grid_fuel_;
  const auto& r_clad = heat_driver_->r_grid_clad_;
  const auto& z = heat_driver_->z_;

  std::unordered_map<CellInstance, int> tracked;

  // index for rings in solid
  int ring_index = 0;

  // TODO: Don't hardcode number of azimuthal segments
  int n_azimuthal = 4;

  int n_fissionable_mapped = 0;

  // Establish mappings between solid regions and OpenMC cells. The center
  // coordinate for each region in the T/H model is obtained and used to
  // determine the OpenMC cell at that position.
  for (gsl::index i = 0; i < heat_driver_->n_pins_; ++i) {
    double x_center = heat_driver_->pin_centers_(i, 0);
    double y_center = heat_driver_->pin_centers_(i, 1);

    for (gsl::index j = 0; j < heat_driver_->n_axial_; ++j) {
      double zavg = 0.5 * (z(j) + z(j + 1));

      for (gsl::index k = 0; k < heat_driver_->n_rings(); ++k) {
        double ravg;
        if (k < heat_driver_->n_fuel_rings_) {
          ravg = 0.5 * (r_fuel(k) + r_fuel(k + 1));
        } else {
          int m = k - heat_driver_->n_fuel_rings_;
          ravg = 0.5 * (r_clad(m) + r_clad(m + 1));
        }

        for (gsl::index m = 0; m < n_azimuthal; ++m) {
          double m_avg = m + 0.5;
          double theta = 2.0 * m_avg * PI / n_azimuthal;
          double x = x_center + ravg * std::cos(theta);
          double y = y_center + ravg * std::sin(theta);

          // Determine cell instance corresponding to given pin location
          Position r{x, y, zavg};
          CellInstance c{r};
          if (tracked.find(c) == tracked.end()) {
            openmc_driver_->cells_.push_back(c);
            tracked[c] = openmc_driver_->cells_.size() - 1;

            if (k < heat_driver_->n_fuel_rings_)
              n_fissionable_mapped++;
          }

          // ensure that the cell being mapped for the pellet region contains
          // a fissionable material as a way to check that the T/H geometry of the
          // pellet region is sufficiently refined to account for all OpenMC cells.
          // This check does not ensure that we have accounted for
          // _all_ fissionable cells, just that the models line up in the pellet region,
          // from which we can infer that the general geometry lines up (because otherwise
          // we could be mapping fluid cells to the surrogate pins, etc.)
          if (k < heat_driver_->n_fuel_rings_) {
            Ensures(c.is_fissionable());
          } else {
            Ensures(!c.is_fissionable());
          }

          // Map OpenMC material to ring and vice versa
          int32_t array_index = tracked[c];
          cell_inst_to_ring_[array_index].push_back(ring_index);
          ring_to_cell_inst_[ring_index].push_back(array_index);
        }

        ++ring_index;
      }
    }
  }

  // check that all fissionable cells have been mapped to the T/H model.
  // Note that this does not check for distortions in the model, just that
  // the surrogate T/H model is sufficiently fine in the pellet region to
  // capture all fissionable cells in OpenMC.
  Ensures(openmc_driver_->n_fissionable_cells_ == n_fissionable_mapped);

  // set the number of solid cells before defining mappings for fluid cells
  if (openmc_driver_->active()) {
    n_solid_cells_ = gsl::narrow<int>(openmc_driver_->cells_.size());
  }

  // index for elements in fluid
  int fluid_index = 0;

  // Establish mappings between fluid regions and OpenMC cells; it is assumed that there
  // are no azimuthal divisions in the fluid phase in the OpenMC model so that we
  // can take a point on a 45 degree ray from the pin center. TODO: add a check to make
  // sure that the T/H model is finer than the OpenMC model.
  int n_elems = heat_driver_->n_pins_ * heat_driver_->n_axial_;
  elem_to_cell_inst_.resize(n_elems);

  for (gsl::index i = 0; i < heat_driver_->n_pins_; ++i) {
    double x_center = heat_driver_->pin_centers_(i, 0);
    double y_center = heat_driver_->pin_centers_(i, 1);

    for (gsl::index j = 0; j < heat_driver_->n_axial_; ++j) {
      double zavg = 0.5 * (z(j) + z(j + 1));
      double l = heat_driver_->pin_pitch() / std::sqrt(2.0);
      double d = (l - heat_driver_->clad_outer_radius_) / 2.0;
      double x = x_center + (heat_driver_->clad_outer_radius_ + d) * std::sqrt(2.0) / 2.0;
      double y = y_center + (heat_driver_->clad_outer_radius_ + d) * std::sqrt(2.0) / 2.0;

      // Determine cell instance corresponding to given fluid location
      Position r{x, y, zavg};
      CellInstance c{r};
      if (tracked.find(c) == tracked.end()) {
        openmc_driver_->cells_.push_back(c);
        tracked[c] = gsl::narrow<int>(openmc_driver_->cells_.size() - 1);

        Ensures(!c.is_fissionable());
      }

      // Map OpenMC material to fluid element and vice versa
      int32_t array_index = tracked[c];
      cell_inst_to_elem_[array_index].push_back(fluid_index);
      elem_to_cell_inst_[fluid_index].push_back(array_index);

      ++fluid_index;
    }
  }

  if (openmc_driver_->active()) {
    n_fluid_cells_ = gsl::narrow<int>(openmc_driver_->cells_.size() - n_solid_cells_);
  }
}

void OpenmcHeatDriver::init_tallies()
{
  using gsl::index;
  using gsl::narrow_cast;

  if (openmc_driver_->active()) {
    // Build vector of cell instances to construct tallies; tallies are only
    // used in the solid regions
    std::vector<openmc::CellInstance> cells;
    for (gsl::index c = 0; c < n_solid_cells_; ++c) {
      const auto& cell = openmc_driver_->cells_[c];
      cells.push_back(
        {narrow_cast<index>(cell.index_), narrow_cast<index>(cell.instance_)});
    }

    openmc_driver_->create_tallies(cells);
  }
}

void OpenmcHeatDriver::init_densities()
{
  if (this->has_global_coupling_data()) {
    std::size_t n_fluid = heat_driver_->n_pins_ * heat_driver_->n_axial_;
    densities_.resize({n_fluid});
    densities_prev_.resize({n_fluid});

    if (density_ic_ == Initial::neutronics) {
      // Loop over all of the fluid regions in the heat transfer model and set
      // the density IC based on the densities used in the OpenMC input file.

      int fluid_index = 0;
      for (gsl::index fluid_index = 0; fluid_index < n_fluid; ++fluid_index) {
        double mass = 0.0;
        double total_vol = 0.0;

        for (const auto& idx : elem_to_cell_inst_[fluid_index]) {
          const auto& c = openmc_driver_->cells_[idx];
          double vol = c.volume_;

          total_vol += vol;
          mass += c.get_density() * vol;
        }

        densities_(fluid_index) = mass / total_vol;
      }

      std::copy(densities_.begin(), densities_.end(), densities_prev_.begin());
    } else if (density_ic_ == Initial::heat) {
      throw std::runtime_error{"Density initial conditions from surrogate heat-fluids "
                               "solver not supported."};
    }
  }
}

void OpenmcHeatDriver::init_temperatures()
{
  if (this->has_global_coupling_data()) {
    std::size_t n_solid =
      heat_driver_->n_pins_ * heat_driver_->n_axial_ * heat_driver_->n_rings();
    std::size_t n_fluid = heat_driver_->n_pins_ * heat_driver_->n_axial_;
    temperatures_.resize({n_solid + n_fluid});
    temperatures_prev_.resize({n_solid + n_fluid});

    if (temperature_ic_ == Initial::neutronics) {
      // Loop over all of the rings in the heat transfer model and set the temperature IC
      // based on temperatures used in the OpenMC input file. More than one OpenMC cell
      // may correspond to a particular ring, so the initial temperature set for that ring
      // should be a volume average of the OpenMC cell temperatures. Then, loop over all
      // axial fluid cells for each pin and set the fluid temperature IC based on the
      // fluid temperature in the OpenMC model. No more than one OpenMC cell should
      // correspond to a given axial fluid cell, so no averaging is needed (though we
      // retain the averaging syntax below: TODO - update elem_to_cell_inst_ to permit
      // multiple OpenMC cells per axial.

      // TODO: This initial condition used in the coupled driver does not truly represent
      // the actual initial condition used in the OpenMC input file, since the surrogate
      // heat solver only includes a radial dependence, though the various azimuthal
      // segments corresponding to a single heat transfer ring may run with different
      // initial temperatures. For this reading of the initial condition to be completely
      // accurate, the heat solver must either be modified to include an azimuthal
      // dependence, or a check inserted here to ensure that the initial temperatures in
      // the azimuthal segments in the OpenMC model are all the same (i.e. no initial
      // azimuthal dependence).

      // set temperatures for solid cells
      for (gsl::index i = 0; i < n_solid; ++i) {
        const auto& cell_instances = ring_to_cell_inst_[i];

        double T_avg = 0.0;
        double total_vol = 0.0;

        for (auto idx : cell_instances) {
          const auto& c = openmc_driver_->cells_[idx];
          double vol = c.volume_;

          total_vol += vol;
          T_avg += c.get_temperature() * vol;
        }

        temperatures_[i] = T_avg / total_vol;
      }

      // set temperatures for fluid cells
      for (gsl::index i = 0; i < n_fluid; ++i) {
        const auto& cell_instances = elem_to_cell_inst_[i];

        double T_avg = 0.0;
        double total_vol = 0.0;

        for (auto idx : cell_instances) {
          const auto& c = openmc_driver_->cells_[idx];
          double vol = c.volume_;

          total_vol += vol;
          T_avg += c.get_temperature() * vol;
        }

        temperatures_[i + n_solid] = T_avg / total_vol;
      }

      std::copy(temperatures_.begin(), temperatures_.end(), temperatures_prev_.begin());
    } else if (temperature_ic_ == Initial::heat) {
      throw std::runtime_error{
        "Temperature initial conditions from surrogate heat-fluids "
        "solver not supported."};
    }
  }
}

void OpenmcHeatDriver::init_heat_source()
{
  heat_source_ = xt::empty<double>({n_solid_cells_});
  heat_source_prev_ = xt::empty<double>({n_solid_cells_});
}

void OpenmcHeatDriver::set_heat_source()
{
  // zero out heat source in single-physics heat solver
  for (auto& val : heat_driver_->source_)
    val = 0.0;

  int ring_index = 0;
  for (gsl::index i = 0; i < heat_driver_->n_pins_; ++i) {
    for (gsl::index j = 0; j < heat_driver_->n_axial_; ++j) {
      // Loop over radial rings
      for (gsl::index k = 0; k < heat_driver_->n_rings(); ++k) {
        // Only update heat source in fuel
        if (k < heat_driver_->n_fuel_rings_) {
          // Get average Q value across each azimuthal segment
          const auto& cell_instances = ring_to_cell_inst_[ring_index];
          double q_avg = 0.0;
          for (auto idx : cell_instances) {
            q_avg += heat_source_(idx);
          }
          q_avg /= cell_instances.size();

          // Set Q in appropriate (pin, axial, ring)
          heat_driver_->source_.at(i, j, k) = q_avg;
        }
        ++ring_index;
      }
    }
  }
}

void OpenmcHeatDriver::set_density()
{
  for (gsl::index i = n_solid_cells_; i < n_solid_cells_ + n_fluid_cells_; ++i) {
    const auto& c = openmc_driver_->cells_[i];
    const auto& elems = cell_inst_to_elem_[i];

    double mass = 0.0;
    double total_vol = 0.0;
    for (auto fluid_index : elems) {
      auto z_index = fluid_index % heat_driver_->n_axial_;
      double vol = heat_driver_->z_(z_index + 1) - heat_driver_->z_(z_index);
      mass += densities_(fluid_index) * vol;
      total_vol += vol;
    }

    c.set_density(mass / total_vol);
  }
}

void OpenmcHeatDriver::set_temperature()
{
  const auto& r_fuel = heat_driver_->r_grid_fuel_;
  const auto& r_clad = heat_driver_->r_grid_clad_;

  // For each OpenMC solid cell, volume average solid temperatures and set
  // by grabbing the T/H regions corresponding to the cell
  for (gsl::index i = 0; i < n_solid_cells_; ++i) {
    const auto& c = openmc_driver_->cells_[i];
    const auto& rings = cell_inst_to_ring_[i];

    double average_temp = 0.0;
    double total_vol = 0.0;
    for (auto solid_index : rings) {
      auto ring_index = solid_index % heat_driver_->n_rings();
      double vol = heat_driver_->solid_areas_(ring_index);
      average_temp += temperatures_(solid_index) * vol;
      total_vol += vol;
    }

    average_temp /= total_vol;
    c.set_temperature(average_temp);
  }

  // For each OpenMC fluid cell, set the fluid temperature by grabbing
  // the T/H region corresponding to the cell
  auto fluid_offset =
    heat_driver_->n_pins_ * heat_driver_->n_rings() * heat_driver_->n_axial_;
  for (gsl::index i = 0; i < n_fluid_cells_; ++i) {
    int index = n_solid_cells_ + i;
    const auto& c = openmc_driver_->cells_[index];
    const auto& elems = cell_inst_to_elem_[index];

    double average_temp = 0.0;
    double total_vol = 0.0;
    for (auto fluid_index : elems) {
      auto z_index = fluid_index % heat_driver_->n_axial_;
      double vol = heat_driver_->z_(z_index + 1) - heat_driver_->z_(z_index);
      average_temp += temperatures_(fluid_index + fluid_offset) * vol;
      total_vol += vol;
    }

    c.set_temperature(average_temp / total_vol);
  }
}

} // namespace enrico
