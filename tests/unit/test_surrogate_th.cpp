/**
 * \file test_surrogate_th.cpp
 * \brief Unit tests for surrogate T/H solver.
 */

#include "catch.hpp"
#include "pugixml.hpp"
#include "enrico/surrogate_heat_driver.h"

TEST_CASE("Verify construction of surrogate thermal-hydraulics driver - single assembly", "[construction]") {
  // load input file
  pugi::xml_document doc;
  auto result = doc.load_file("inputs/test_surrogate_th_single.xml");

  CHECK(result);

  auto root = doc.document_element();
  auto node = root.child("heat_fluids");

  enrico::SurrogateHeatDriver driver(MPI_COMM_NULL, node);
  enrico::SurrogateHeatDriverAssembly assembly = driver.assembly_drivers_[0];

  SECTION("Verify XML input file reading") {
    CHECK(driver.clad_inner_radius() == Approx(0.414));
    CHECK(driver.clad_outer_radius() == Approx(0.475));
    CHECK(driver.pellet_radius() == Approx(0.406));
    CHECK(driver.n_fuel_rings() == 5);
    CHECK(driver.n_clad_rings() == 3);
    CHECK(driver.n_pins_x() == 7);
    CHECK(driver.n_pins_y() == 4);
    CHECK(driver.pin_pitch() == Approx(1.26));
    CHECK(driver.inlet_temperature() == Approx(500.0));
    CHECK(driver.mass_flowrate() == Approx(0.3));
  }

  SECTION("Verify setting of default XML inputs") {
    CHECK(driver.max_subchannel_its() == 100);
    CHECK(driver.subchannel_tol_h() == Approx(1.0e-2));
    CHECK(driver.subchannel_tol_p() == Approx(1.0e-2));
    CHECK(driver.heat_tol() == Approx(1.0e-4));
  }

  SECTION("Verify calculation of pin center coordinates") {
    const auto& centers = assembly.pin_centers_;

    // check y-coordinates row by row
    for (int i = 0; i < 7; ++i)
      CHECK(assembly.pin_centers_(i, 1) == Approx(1.89));
    for (int i = 7; i < 14; ++i)
      CHECK(centers(i, 1) == Approx(0.63));
    for (int i = 14; i < 21; ++i)
      CHECK(centers(i, 1) == Approx(-0.63));
    for (int i = 21; i < 28; ++i)
      CHECK(centers(i, 1) == Approx(-1.89));

    // check x-coordinates column by column
    for (int i = 0; i < 22; i += 7)
      CHECK(centers(i, 0) == Approx(-3.78));
    for (int i = 1; i < 23; i += 7)
      CHECK(centers(i, 0) == Approx(-2.52));
    for (int i = 2; i < 24; i += 7)
      CHECK(centers(i, 0) == Approx(-1.26));
    for (int i = 3; i < 25; i += 7)
      CHECK(centers(i, 0) == Approx(0.0));
    for (int i = 4; i < 26; i += 7)
      CHECK(centers(i, 0) == Approx(1.26));
    for (int i = 5; i < 27; i += 7)
      CHECK(centers(i, 0) == Approx(2.52));
    for (int i = 6; i < 28; i += 7)
      CHECK(centers(i, 0) == Approx(3.78));
  }

  SECTION("Verify construction of channels") {
    const auto& channels = assembly.channels_;

    CHECK(channels.size() == 40);

    CHECK(channels[0].rod_ids_[0] == 0);

    CHECK(channels[1].rod_ids_[0] == 0);
    CHECK(channels[1].rod_ids_[1] == 1);

    CHECK(channels[2].rod_ids_[0] == 1);
    CHECK(channels[2].rod_ids_[1] == 2);

    CHECK(channels[3].rod_ids_[0] == 2);
    CHECK(channels[3].rod_ids_[1] == 3);

    CHECK(channels[4].rod_ids_[0] == 3);
    CHECK(channels[4].rod_ids_[1] == 4);

    CHECK(channels[5].rod_ids_[0] == 4);
    CHECK(channels[5].rod_ids_[1] == 5);

    CHECK(channels[6].rod_ids_[0] == 5);
    CHECK(channels[6].rod_ids_[1] == 6);

    CHECK(channels[7].rod_ids_[0] == 6);

    CHECK(channels[8].rod_ids_[0] == 0);
    CHECK(channels[8].rod_ids_[1] == 7);

    CHECK(channels[9].rod_ids_[0] == 0);
    CHECK(channels[9].rod_ids_[1] == 1);
    CHECK(channels[9].rod_ids_[2] == 7);
    CHECK(channels[9].rod_ids_[3] == 8);

    CHECK(channels[10].rod_ids_[0] == 1);
    CHECK(channels[10].rod_ids_[1] == 2);
    CHECK(channels[10].rod_ids_[2] == 8);
    CHECK(channels[10].rod_ids_[3] == 9);

    CHECK(channels[11].rod_ids_[0] == 2);
    CHECK(channels[11].rod_ids_[1] == 3);
    CHECK(channels[11].rod_ids_[2] == 9);
    CHECK(channels[11].rod_ids_[3] == 10);

    CHECK(channels[12].rod_ids_[0] == 3);
    CHECK(channels[12].rod_ids_[1] == 4);
    CHECK(channels[12].rod_ids_[2] == 10);
    CHECK(channels[12].rod_ids_[3] == 11);

    CHECK(channels[13].rod_ids_[0] == 4);
    CHECK(channels[13].rod_ids_[1] == 5);
    CHECK(channels[13].rod_ids_[2] == 11);
    CHECK(channels[13].rod_ids_[3] == 12);

    CHECK(channels[14].rod_ids_[0] == 5);
    CHECK(channels[14].rod_ids_[1] == 6);
    CHECK(channels[14].rod_ids_[2] == 12);
    CHECK(channels[14].rod_ids_[3] == 13);

    CHECK(channels[15].rod_ids_[0] == 6);
    CHECK(channels[15].rod_ids_[1] == 13);

    CHECK(channels[16].rod_ids_[0] == 7);
    CHECK(channels[16].rod_ids_[1] == 14);

    CHECK(channels[17].rod_ids_[0] == 7);
    CHECK(channels[17].rod_ids_[1] == 8);
    CHECK(channels[17].rod_ids_[2] == 14);
    CHECK(channels[17].rod_ids_[3] == 15);

    CHECK(channels[18].rod_ids_[0] == 8);
    CHECK(channels[18].rod_ids_[1] == 9);
    CHECK(channels[18].rod_ids_[2] == 15);
    CHECK(channels[18].rod_ids_[3] == 16);

    CHECK(channels[19].rod_ids_[0] == 9);
    CHECK(channels[19].rod_ids_[1] == 10);
    CHECK(channels[19].rod_ids_[2] == 16);
    CHECK(channels[19].rod_ids_[3] == 17);

    CHECK(channels[20].rod_ids_[0] == 10);
    CHECK(channels[20].rod_ids_[1] == 11);
    CHECK(channels[20].rod_ids_[2] == 17);
    CHECK(channels[20].rod_ids_[3] == 18);

    CHECK(channels[21].rod_ids_[0] == 11);
    CHECK(channels[21].rod_ids_[1] == 12);
    CHECK(channels[21].rod_ids_[2] == 18);
    CHECK(channels[21].rod_ids_[3] == 19);

    CHECK(channels[22].rod_ids_[0] == 12);
    CHECK(channels[22].rod_ids_[1] == 13);
    CHECK(channels[22].rod_ids_[2] == 19);
    CHECK(channels[22].rod_ids_[3] == 20);

    CHECK(channels[23].rod_ids_[0] == 13);
    CHECK(channels[23].rod_ids_[1] == 20);

    CHECK(channels[24].rod_ids_[0] == 14);
    CHECK(channels[24].rod_ids_[1] == 21);

    CHECK(channels[25].rod_ids_[0] == 14);
    CHECK(channels[25].rod_ids_[1] == 15);
    CHECK(channels[25].rod_ids_[2] == 21);
    CHECK(channels[25].rod_ids_[3] == 22);

    CHECK(channels[26].rod_ids_[0] == 15);
    CHECK(channels[26].rod_ids_[1] == 16);
    CHECK(channels[26].rod_ids_[2] == 22);
    CHECK(channels[26].rod_ids_[3] == 23);

    CHECK(channels[27].rod_ids_[0] == 16);
    CHECK(channels[27].rod_ids_[1] == 17);
    CHECK(channels[27].rod_ids_[2] == 23);
    CHECK(channels[27].rod_ids_[3] == 24);

    CHECK(channels[28].rod_ids_[0] == 17);
    CHECK(channels[28].rod_ids_[1] == 18);
    CHECK(channels[28].rod_ids_[2] == 24);
    CHECK(channels[28].rod_ids_[3] == 25);

    CHECK(channels[29].rod_ids_[0] == 18);
    CHECK(channels[29].rod_ids_[1] == 19);
    CHECK(channels[29].rod_ids_[2] == 25);
    CHECK(channels[29].rod_ids_[3] == 26);

    CHECK(channels[30].rod_ids_[0] == 19);
    CHECK(channels[30].rod_ids_[1] == 20);
    CHECK(channels[30].rod_ids_[2] == 26);
    CHECK(channels[30].rod_ids_[3] == 27);

    CHECK(channels[31].rod_ids_[0] == 20);
    CHECK(channels[31].rod_ids_[1] == 27);

    CHECK(channels[32].rod_ids_[0] == 21);

    CHECK(channels[33].rod_ids_[0] == 21);
    CHECK(channels[33].rod_ids_[1] == 22);

    CHECK(channels[34].rod_ids_[0] == 22);
    CHECK(channels[34].rod_ids_[1] == 23);

    CHECK(channels[35].rod_ids_[0] == 23);
    CHECK(channels[35].rod_ids_[1] == 24);

    CHECK(channels[36].rod_ids_[0] == 24);
    CHECK(channels[36].rod_ids_[1] == 25);

    CHECK(channels[37].rod_ids_[0] == 25);
    CHECK(channels[37].rod_ids_[1] == 26);

    CHECK(channels[38].rod_ids_[0] == 26);
    CHECK(channels[38].rod_ids_[1] == 27);

    CHECK(channels[39].rod_ids_[0] == 27);

    std::vector<int> corner = {0, 7, 32, 39};
    std::vector<int> edge = {1, 2, 3, 4, 5, 6, 8, 15, 16, 23, 24, 31, 33, 34, 35, 36, 37, 38};
    double corner_area = 0.21969453938345077;
    double interior_area = 0.8787781575338031;
    double edge_area = 0.43938907876;
    double total_area = 4.0 * corner_area + 18.0 * interior_area + 18.0 * edge_area;

    for (int i = 0; i < 40; ++i) {
      if (std::count(corner.begin(), corner.end(), i)) {
        CHECK(channels[i].area_ == Approx(corner_area));
        CHECK(assembly.channel_flowrates_(i) == Approx(corner_area / total_area * 0.3));
      }
      else {
        if (std::count(edge.begin(), edge.end(), i)) {
          CHECK(channels[i].area_ == Approx(edge_area));
          CHECK(assembly.channel_flowrates_(i) == Approx(edge_area / total_area * 0.3));
        }
        else {
          CHECK(channels[i].area_ == Approx(interior_area));
          CHECK(assembly.channel_flowrates_(i) == Approx(interior_area / total_area * 0.3));
        }
      }
    }
  }

  SECTION("Verify construction of rods") {
    const auto& rods = assembly.rods_;

    CHECK(rods.size() == 28);

    for (int i = 0; i < rods.size(); ++i) {
      CHECK(rods[i].clad_outer_radius_ == Approx(0.475));
      CHECK(rods[i].clad_inner_radius_ == Approx(0.414));
      CHECK(rods[i].pellet_radius_ == Approx(0.406));
      CHECK(rods[i].channel_ids_.size() == 4);
    }

    CHECK(rods[0].channel_ids_[0] == 0);
    CHECK(rods[0].channel_ids_[1] == 1);
    CHECK(rods[0].channel_ids_[2] == 8);
    CHECK(rods[0].channel_ids_[3] == 9);

    CHECK(rods[1].channel_ids_[0] == 1);
    CHECK(rods[1].channel_ids_[1] == 2);
    CHECK(rods[1].channel_ids_[2] == 9);
    CHECK(rods[1].channel_ids_[3] == 10);

    CHECK(rods[2].channel_ids_[0] == 2);
    CHECK(rods[2].channel_ids_[1] == 3);
    CHECK(rods[2].channel_ids_[2] == 10);
    CHECK(rods[2].channel_ids_[3] == 11);

    CHECK(rods[3].channel_ids_[0] == 3);
    CHECK(rods[3].channel_ids_[1] == 4);
    CHECK(rods[3].channel_ids_[2] == 11);
    CHECK(rods[3].channel_ids_[3] == 12);

    CHECK(rods[4].channel_ids_[0] == 4);
    CHECK(rods[4].channel_ids_[1] == 5);
    CHECK(rods[4].channel_ids_[2] == 12);
    CHECK(rods[4].channel_ids_[3] == 13);

    CHECK(rods[5].channel_ids_[0] == 5);
    CHECK(rods[5].channel_ids_[1] == 6);
    CHECK(rods[5].channel_ids_[2] == 13);
    CHECK(rods[5].channel_ids_[3] == 14);

    CHECK(rods[6].channel_ids_[0] == 6);
    CHECK(rods[6].channel_ids_[1] == 7);
    CHECK(rods[6].channel_ids_[2] == 14);
    CHECK(rods[6].channel_ids_[3] == 15);

    CHECK(rods[7].channel_ids_[0] == 8);
    CHECK(rods[7].channel_ids_[1] == 9);
    CHECK(rods[7].channel_ids_[2] == 16);
    CHECK(rods[7].channel_ids_[3] == 17);

    CHECK(rods[8].channel_ids_[0] == 9);
    CHECK(rods[8].channel_ids_[1] == 10);
    CHECK(rods[8].channel_ids_[2] == 17);
    CHECK(rods[8].channel_ids_[3] == 18);

    CHECK(rods[9].channel_ids_[0] == 10);
    CHECK(rods[9].channel_ids_[1] == 11);
    CHECK(rods[9].channel_ids_[2] == 18);
    CHECK(rods[9].channel_ids_[3] == 19);

    CHECK(rods[10].channel_ids_[0] == 11);
    CHECK(rods[10].channel_ids_[1] == 12);
    CHECK(rods[10].channel_ids_[2] == 19);
    CHECK(rods[10].channel_ids_[3] == 20);

    CHECK(rods[11].channel_ids_[0] == 12);
    CHECK(rods[11].channel_ids_[1] == 13);
    CHECK(rods[11].channel_ids_[2] == 20);
    CHECK(rods[11].channel_ids_[3] == 21);

    CHECK(rods[12].channel_ids_[0] == 13);
    CHECK(rods[12].channel_ids_[1] == 14);
    CHECK(rods[12].channel_ids_[2] == 21);
    CHECK(rods[12].channel_ids_[3] == 22);

    CHECK(rods[13].channel_ids_[0] == 14);
    CHECK(rods[13].channel_ids_[1] == 15);
    CHECK(rods[13].channel_ids_[2] == 22);
    CHECK(rods[13].channel_ids_[3] == 23);

    CHECK(rods[14].channel_ids_[0] == 16);
    CHECK(rods[14].channel_ids_[1] == 17);
    CHECK(rods[14].channel_ids_[2] == 24);
    CHECK(rods[14].channel_ids_[3] == 25);

    CHECK(rods[15].channel_ids_[0] == 17);
    CHECK(rods[15].channel_ids_[1] == 18);
    CHECK(rods[15].channel_ids_[2] == 25);
    CHECK(rods[15].channel_ids_[3] == 26);

    CHECK(rods[16].channel_ids_[0] == 18);
    CHECK(rods[16].channel_ids_[1] == 19);
    CHECK(rods[16].channel_ids_[2] == 26);
    CHECK(rods[16].channel_ids_[3] == 27);

    CHECK(rods[17].channel_ids_[0] == 19);
    CHECK(rods[17].channel_ids_[1] == 20);
    CHECK(rods[17].channel_ids_[2] == 27);
    CHECK(rods[17].channel_ids_[3] == 28);

    CHECK(rods[18].channel_ids_[0] == 20);
    CHECK(rods[18].channel_ids_[1] == 21);
    CHECK(rods[18].channel_ids_[2] == 28);
    CHECK(rods[18].channel_ids_[3] == 29);

    CHECK(rods[19].channel_ids_[0] == 21);
    CHECK(rods[19].channel_ids_[1] == 22);
    CHECK(rods[19].channel_ids_[2] == 29);
    CHECK(rods[19].channel_ids_[3] == 30);

    CHECK(rods[20].channel_ids_[0] == 22);
    CHECK(rods[20].channel_ids_[1] == 23);
    CHECK(rods[20].channel_ids_[2] == 30);
    CHECK(rods[20].channel_ids_[3] == 31);

    CHECK(rods[21].channel_ids_[0] == 24);
    CHECK(rods[21].channel_ids_[1] == 25);
    CHECK(rods[21].channel_ids_[2] == 32);
    CHECK(rods[21].channel_ids_[3] == 33);

    CHECK(rods[22].channel_ids_[0] == 25);
    CHECK(rods[22].channel_ids_[1] == 26);
    CHECK(rods[22].channel_ids_[2] == 33);
    CHECK(rods[22].channel_ids_[3] == 34);

    CHECK(rods[23].channel_ids_[0] == 26);
    CHECK(rods[23].channel_ids_[1] == 27);
    CHECK(rods[23].channel_ids_[2] == 34);
    CHECK(rods[23].channel_ids_[3] == 35);

    CHECK(rods[24].channel_ids_[0] == 27);
    CHECK(rods[24].channel_ids_[1] == 28);
    CHECK(rods[24].channel_ids_[2] == 35);
    CHECK(rods[24].channel_ids_[3] == 36);

    CHECK(rods[25].channel_ids_[0] == 28);
    CHECK(rods[25].channel_ids_[1] == 29);
    CHECK(rods[25].channel_ids_[2] == 36);
    CHECK(rods[25].channel_ids_[3] == 37);

    CHECK(rods[26].channel_ids_[0] == 29);
    CHECK(rods[26].channel_ids_[1] == 30);
    CHECK(rods[26].channel_ids_[2] == 37);
    CHECK(rods[26].channel_ids_[3] == 38);

    CHECK(rods[27].channel_ids_[0] == 30);
    CHECK(rods[27].channel_ids_[1] == 31);
    CHECK(rods[27].channel_ids_[2] == 38);
    CHECK(rods[27].channel_ids_[3] == 39);
  }

  SECTION("Verify axial discretization") {
    const auto& z = driver.z_;

    CHECK(driver.n_axial_ == 6);

    CHECK(z(0) == Approx(0.1));
    CHECK(z(1) == Approx(0.5));
    CHECK(z(2) == Approx(1.1));
    CHECK(z(3) == Approx(1.4));
    CHECK(z(4) == Approx(2.0));
    CHECK(z(5) == Approx(2.1));
    CHECK(z(6) == Approx(2.2));
  }

  SECTION("Verify radial discretization") {
    const auto& rf = assembly.r_grid_fuel_;

    CHECK(rf(0) == Approx(0.0));
    CHECK(rf(1) == Approx(0.0812));
    CHECK(rf(2) == Approx(0.1624));
    CHECK(rf(3) == Approx(0.2436));
    CHECK(rf(4) == Approx(0.3248));
    CHECK(rf(5) == Approx(0.406));

    const auto& rc = assembly.r_grid_clad_;

    CHECK(rc(0) == Approx(0.414));
    CHECK(rc(1) == Approx(0.434333333333333));
    CHECK(rc(2) == Approx(0.454666666666666));
    CHECK(rc(3) == Approx(0.475));
  }

}
