/**
 * @author Julius Isken, Max Engel
 *
 * @section LICENSE
 * Copyright 2020, Julius Isken
 *
 * Redistribution and use in source and binary forms, with or without
 *modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *this list of conditions and the following disclaimer in the documentation
 *and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Tests the Shock-Shock setup.
 **/
#include "ShockShock1d.h"
#include <catch2/catch.hpp>

TEST_CASE("Test the 1-dimensional shock shock setup", "[ShockShock1d]") {
  tsunami_lab::setups::ShockShock1d l_shockShock(10, 15, 20);

  // Check the momentum and the height of the setup
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
