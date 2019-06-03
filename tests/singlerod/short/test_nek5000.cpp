#include "enrico/nek_driver.h"
#include <mpi.h>

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  // Parse enrico.xml file
  pugi::xml_document doc;
  auto result = doc.load_file("enrico.xml");
  if (!result) {
    throw std::runtime_error{"Unable to load enrico.xml file"};
  }

  // Get root element
  auto root = doc.document_element();

  {
    enrico::NekDriver test_driver(MPI_COMM_WORLD,
                                  root.child("pressure_bc").text().as_double(),
                                  root.child("nek5000"));
    test_driver.init_step();
    test_driver.solve_step();
    test_driver.finalize_step();
  }

  MPI_Finalize();

  return 0;
}
