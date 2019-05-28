#include "enrico/openmc_heat_driver.h"

#include "enrico/error.h"
#include "enrico/message_passing.h"
#include "openmc/constants.h"
#include "xtensor/xstrided_view.hpp"

#include <cmath>
#include <gsl/gsl>
#include <unordered_map>

namespace enrico {

OpenmcHeatDriver::OpenmcHeatDriver(MPI_Comm comm, pugi::xml_node node)
  : CoupledDriver(comm, node)
{
  // Initialize OpenMC and surrogate heat drivers
  openmc_driver_ = std::make_unique<OpenmcDriver>(comm);
  pugi::xml_node surr_node = node.child("heat_surrogate");
  heat_driver_ = std::make_unique<SurrogateHeatDriver>(comm, surr_node);

  // Create mappings for fuel pins and setup tallies for OpenMC
  init_mappings();
  init_tallies();

  init_temperatures();
  init_heat_source();
}

NeutronicsDriver& OpenmcHeatDriver::get_neutronics_driver() const
{
  return *openmc_driver_;
}

Driver& OpenmcHeatDriver::get_heat_driver() const
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
  int ring_index = 0;
  // TODO: Don't hardcode number of azimuthal segments
  int n_azimuthal = 4;

  for (int i = 0; i < heat_driver_->n_pins_; ++i) {
    for (int j = 0; j < heat_driver_->n_axial_; ++j) {
      // Get average z value
      double zavg = 0.5 * (z(j) + z(j + 1));

      // Loop over radial rings
      for (int k = 0; k < heat_driver_->n_rings(); ++k) {
        double ravg;
        if (k < heat_driver_->n_fuel_rings_) {
          ravg = 0.5 * (r_fuel(k) + r_fuel(k + 1));
        } else {
          int m = k - heat_driver_->n_fuel_rings_;
          ravg = 0.5 * (r_clad(m) + r_clad(m + 1));
        }

        for (int m = 0; m < n_azimuthal; ++m) {
          double m_avg = m + 0.5;
          double theta = 2.0 * m_avg * PI / n_azimuthal;
          double x = ravg * std::cos(theta);
          double y = ravg * std::sin(theta);

          // Determine cell instance corresponding to given pin location
          Position r{x, y, zavg};
          CellInstance c{r};
          if (tracked.find(c) == tracked.end()) {
            openmc_driver_->cells_.push_back(c);
            tracked[c] = openmc_driver_->cells_.size() - 1;
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

  if (openmc_driver_->active()) {
    n_materials_ = openmc_driver_->cells_.size();
  }
}

void OpenmcHeatDriver::init_tallies()
{
  if (openmc_driver_->active()) {
    // Build vector of material indices
    std::vector<int32_t> mats;
    for (const auto& c : openmc_driver_->cells_) {
      mats.push_back(c.material_index_);
    }
    openmc_driver_->create_tallies(mats);
  }
}

void OpenmcHeatDriver::init_temperatures()
{
  std::size_t n =
    heat_driver_->n_pins_ * heat_driver_->n_axial_ * heat_driver_->n_rings();
  temperatures_.resize({n});
  temperatures_prev_.resize({n});

  if (temperature_ic_ == Initial::neutronics) {
    // Loop over all of the rings in the heat transfer model and set the temperature IC
    // based on temperatures used in the OpenMC input file. More than one OpenMC cell may
    // correspond to a particular ring, so the initial temperature set for that ring
    // should be a volume average of the OpenMC cell temperatures.

    // TODO: This initial condition used in the coupled driver does not truly represent
    // the actual initial condition used in the OpenMC input file, since the surrogate
    // heat solver only includes a radial dependence, though the various azimuthal
    // segments corresponding to a single heat transfer ring may run with different
    // initial temperatures. For this reading of the initial condition to be completely
    // accurate, the heat solver must either be modified to include an azimuthal
    // dependence, or a check inserted here to ensure that the initial temperatures in the
    // azimuthal segments in the OpenMC model are all the same (i.e. no initial azimuthal
    // dependence).

    int ring_index = 0;
    for (int i = 0; i < heat_driver_->n_pins_; ++i) {
      for (int j = 0; j < heat_driver_->n_axial_; ++j) {
        for (int k = 0; k < heat_driver_->n_rings(); ++k) {
          const auto& cell_instances = ring_to_cell_inst_[ring_index];

          double T_avg = 0.0;
          double total_vol = 0.0;

          for (auto idx : cell_instances) {
            const auto& c = openmc_driver_->cells_[idx];
            double vol = c.volume_;

            total_vol += vol;
            T_avg += c.get_temperature() * vol;
          }

          T_avg /= total_vol;

          int t_index = (i * heat_driver_->n_axial_ + j) * heat_driver_->n_rings() + k;
          temperatures_prev_[t_index] = T_avg;
          temperatures_[t_index] = T_avg;

          ++ring_index;
        }
      }
    }
  }

  if (temperature_ic_ == Initial::heat) {
    throw std::runtime_error{"Temperature initial conditions from surrogate heat-fluids "
                             "solver not supported."};
  }
}

void OpenmcHeatDriver::init_heat_source()
{
  heat_source_ = xt::empty<double>({n_materials_});
  heat_source_prev_ = xt::empty<double>({n_materials_});
}

void OpenmcHeatDriver::set_heat_source()
{
  // zero out heat source in single-physics heat solver
  for (auto& val : heat_driver_->source_)
    val = 0.0;

  int ring_index = 0;
  for (int i = 0; i < heat_driver_->n_pins_; ++i) {
    for (int j = 0; j < heat_driver_->n_axial_; ++j) {
      // Loop over radial rings
      for (int k = 0; k < heat_driver_->n_rings(); ++k) {
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

void OpenmcHeatDriver::update_temperature()
{
  std::copy(temperatures_.begin(), temperatures_.end(), temperatures_prev_.begin());

  const auto& r_fuel = heat_driver_->r_grid_fuel_;
  const auto& r_clad = heat_driver_->r_grid_clad_;

  temperatures_ = heat_driver_->temperature();

  // For each OpenMC material, volume average temperatures and set
  for (int i = 0; i < openmc_driver_->cells_.size(); ++i) {
    // Get cell instance
    const auto& c = openmc_driver_->cells_[i];

    // Get rings corresponding to this cell instance
    const auto& rings = cell_inst_to_ring_[i];

    // Get volume-average temperature for this material
    double average_temp = 0.0;
    double total_vol = 0.0;
    for (int ring_index : rings) {
      // Use difference in r**2 as a proxy for volume. This is only used for
      // averaging, so the absolute value doesn't matter
      int j = ring_index % heat_driver_->n_rings();
      double vol;
      if (j < heat_driver_->n_fuel_rings_) {
        vol = r_fuel(j + 1) * r_fuel(j + 1) - r_fuel(j) * r_fuel(j);
      } else {
        j -= heat_driver_->n_fuel_rings_;
        vol = r_clad(j + 1) * r_clad(j + 1) - r_clad(j) * r_clad(j);
      }

      average_temp += temperatures_(ring_index) * vol;
      total_vol += vol;
    }
    average_temp /= total_vol;

    // Set temperature for cell instance
    c.set_temperature(average_temp);
  }
}

bool OpenmcHeatDriver::is_converged()
{
  bool converged;
  if (comm_.rank == 0) {
    double norm;
    compute_temperature_norm(Norm::LINF, norm, converged);

    std::string msg = "temperature norm_linf: " + std::to_string(norm);
    comm_.message(msg);
  }
  err_chk(comm_.Bcast(&converged, 1, MPI_CXX_BOOL));
  return converged;
}

} // namespace enrico
