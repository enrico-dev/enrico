#include <stdexcept>

#include "pugixml.hpp"
#include <mpi.h>

#include "enrico/coupled_driver.h"
#include "enrico/mpi_types.h"

int main(int argc, char* argv[])
{
  // Initialize MPI
  MPI_Init(&argc, &argv);
  enrico::init_mpi_datatypes();

  // Define enums for selecting drivers
  enum class Transport { OpenMC, Shift, Surrogate };

  // Parse enrico.xml file
  pugi::xml_document doc;
  auto result = doc.load_file("enrico.xml");
  if (!result) {
    throw std::runtime_error{"Unable to load enrico.xml file"};
  }

  // Get root element
  auto root = doc.document_element();

  // Determine transport driver
  auto neut_driver = std::string{root.child("neutronics").child_value("driver")};
  auto heat_driver = std::string{root.child("heat_fluids").child_value("driver")};
  Transport driver_transport;
  if (neut_driver == "openmc") {
    driver_transport = Transport::OpenMC;
  } else if (neut_driver == "shift") {
    driver_transport = Transport::Shift;
  } else if (neut_driver == "surrogate") {
    driver_transport = Transport::Surrogate;
  } else {
    throw std::runtime_error{"Invalid value for <neutronics><driver>"};
  }

  // Create driver according to selections
  switch (driver_transport) {
  case Transport::OpenMC:
  case Transport::Shift: {
    enrico::CoupledDriver driver{MPI_COMM_WORLD, root};
    driver.execute();

    driver.comm_.message("CoupledDriver times (seconds)");
    driver.timer_report();

    driver.comm_.Barrier();
    driver.comm_.message("NeutronicsDriver times (seconds)");
    driver.get_neutronics_driver().timer_report();

    driver.comm_.Barrier();
    driver.comm_.message("HeatFluidsDriver times (seconds)");
    driver.get_heat_driver().timer_report();

  } break;
  case Transport::Surrogate:
    throw std::runtime_error{"No surrogate particle transport driver implemented"};
    break;
  }

  enrico::free_mpi_datatypes();
  MPI_Finalize();
  return 0;
}
