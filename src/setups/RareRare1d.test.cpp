#include "RareRare1d.h"
#include <catch2/catch.hpp>

TEST_CASE("Test the 1-dimensional rare rare setup", "[RareRare1d]") {
  tsunami_lab::setups::RareRare1d l_rareRare(10, 15, 20);

  // Check epicenter
  REQUIRE(l_rareRare.getHeight(2, 0) == 20);

  REQUIRE(l_rareRare.getMomentumX(2, 0) == -10);

  REQUIRE(l_rareRare.getMomentumY(2, 0) == 0);

  REQUIRE(l_rareRare.getHeight(2, 5) == 20);

  REQUIRE(l_rareRare.getMomentumX(2, 5) == -10);

  REQUIRE(l_rareRare.getMomentumY(2, 2) == 0);

  REQUIRE(l_rareRare.getHeight(20, 0) == 20);

  REQUIRE(l_rareRare.getMomentumX(20, 0) == 10);

  REQUIRE(l_rareRare.getMomentumY(20, 0) == 0);
}
