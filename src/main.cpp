#include <stdexcept>

#include "pugixml.hpp"
#include <mpi.h>

#include "enrico/mpi_types.h"
#include "enrico/openmc_nek_driver.h"

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
  auto s = std::string{root.child_value("driver_transport")};
  Transport driver_transport;
  if (s == "openmc") {
    driver_transport = Transport::OpenMC;
  } else if (s == "surrogate") {
    driver_transport = Transport::Surrogate;
  } else {
    throw std::runtime_error{"Invalid value for <driver_transport>"};
  }

  // Create driver according to selections
  switch (driver_transport) {
  case Transport::OpenMC: {
    enrico::OpenmcNekDriver driver{MPI_COMM_WORLD, root};
    driver.execute();
  } break;
  case Transport::Shift:
    throw std::runtime_error{"Shift transport driver not implemented"};
    break;
  case Transport::Surrogate:
    throw std::runtime_error{"No surrogate particle transport driver implemented"};
    break;
  }

  enrico::free_mpi_datatypes();
  MPI_Finalize();
  return 0;
}
