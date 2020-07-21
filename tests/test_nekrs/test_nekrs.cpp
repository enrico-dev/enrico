#include "enrico/nekrs_driver.h"

int main(int argc, char* argv[]) {
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
    enrico::NekRSDriver test_driver(MPI_COMM_WORLD, root.child("heat_fluids"));
    test_driver.init_step();
    test_driver.solve_step();
    test_driver.finalize_step();
  }

  MPI_Finalize();

  return 0;
}