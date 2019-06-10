#include "smrt/smrt_coupled_driver.h"
#include "pugixml.hpp"
#include <stdexcept>

namespace enrico {

SmrtCoupledDriver::SmrtCoupledDriver()
{
  // TODO: Belongs to main.cpp
  pugi::xml_document doc;
  auto result = doc.load_file("enrico.xml");
  if (!result) {
    throw std::runtime_error{"Unable to load enrico.xml file"};
  }

  // Get root element
  auto root = doc.document_element();

  power_ = root.child("power").text().as_double();
  max_picard_iter_ = root.child("max_picard_iter").text().as_int();
}

} // end namespace enrico
