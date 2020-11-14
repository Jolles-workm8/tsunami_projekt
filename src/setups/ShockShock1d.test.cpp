#include "ShockShock1d.h"
#include <catch2/catch.hpp>

TEST_CASE("Test the 1-dimensional shock shock setup", "[ShockShock1d]") {
  tsunami_lab::setups::ShockShock1d l_shockShock(10, 15, 20);

  // Check epicenter
  REQUIRE(l_shockShock.getHeight(2, 0) == 20);

  REQUIRE(l_shockShock.getMomentumX(2, 0) == 10);

  REQUIRE(l_shockShock.getMomentumY(2, 0) == 0);

  REQUIRE(l_shockShock.getHeight(2, 5) == 20);

  REQUIRE(l_shockShock.getMomentumX(2, 5) == 10);

  REQUIRE(l_shockShock.getMomentumY(2, 2) == 0);

  REQUIRE(l_shockShock.getHeight(20, 0) == 20);

  REQUIRE(l_shockShock.getMomentumX(20, 0) == -10);

  REQUIRE(l_shockShock.getMomentumY(20, 0) == 0);
}
