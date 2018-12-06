#include <stdexcept>

#include <mpi.h>
#include "pugixml.hpp"

#include "stream/message_passing.h"
#include "stream/openmc_nek_driver.h"
#include "stream/openmc_heat_driver.h"

int main(int argc, char* argv[])
{
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Define enums for selecting drivers
  enum class Transport { OpenMC, Shift, Surrogate };
  enum class HeatFluids { Nek5000, Surrogate };

  // Parse stream.xml file
  pugi::xml_document doc;
  auto result = doc.load_file("stream.xml");
  if (!result) {
    throw std::runtime_error{"Unable to load stream.xml file"};
  }

  // Get root element
  auto root = doc.document_element();

  // Determine transport driver
  auto s = std::string{root.child_value("driver_transport")};
  Transport driver_transport;
  if (s == "openmc") {
    driver_transport = Transport::OpenMC;
  } else if (s == "surrogate") {
    driver_transport = Transport::Surrogate;
  } else {
    throw std::runtime_error{"Invalid value for <driver_transport>"};
  }

  // Determine heat/fluids driver
  s = std::string{root.child_value("driver_heatfluids")};
  HeatFluids driver_heatfluids;
  if (s == "nek5000") {
    driver_heatfluids = HeatFluids::Nek5000;
  } else if (s == "surrogate") {
    driver_heatfluids = HeatFluids::Surrogate;
  } else {
    throw std::runtime_error{"Invalid value for <driver_heatfluids>"};
  }

  // Determine power
  double power = root.child("power").text().as_double();

  // Create driver according to selections
  switch (driver_transport) {
  case Transport::OpenMC:
    switch (driver_heatfluids) {
    case HeatFluids::Nek5000:
      {
        stream::OpenmcNekDriver driver {MPI_COMM_WORLD, root};
        driver.solve_in_time();
      }
      break;
    case HeatFluids::Surrogate:
      {
        // Pass XML node for reading settings
        stream::OpenmcHeatDriver driver {MPI_COMM_WORLD, root};
        driver.solve_step();
        driver.to_vtk();
      }
      break;
    }
    break;
  case Transport::Shift:
    throw std::runtime_error{"Shift transport driver not implemented"};
    break;
  case Transport::Surrogate:
    throw std::runtime_error{"No surrogate particle transport driver implemented"};
    break;
  }

  MPI_Finalize();
  return 0;
}
