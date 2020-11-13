#include <catch2/catch.hpp>
#define private public
#include "fwave.h"
#undef public

TEST_CASE("Test Eigenvalues lambda", "[lambdas]") {
  /*
   * Test case:
   *  h: 10 | 9
   *  hu: -30 | 27
   *  u: -3| 3
   *
   * roe height: 9.5
   * roe velocity: (sqrt(10) * -3 + 3 * 3) / ( sqrt(10) + sqrt(9) )
   * = -0.07900211696917202
   * sqrt(g*height): 3.131557121*sqrt(9.5) =  9.652107283894807
   *
   * lambda[0]: -0.07900211696917202 - 9.652107283894807 = -9.731109400863979
   * lambda[1]: -0.07900211696917202 + 9.652107283894807 =  9.573105166925635
   *
   *
   */

  float l_waveSpeedL = 0;
  float l_waveSpeedR = 0;
  tsunami_lab::solvers::fwave::waveSpeeds(10, 9, -3, 3, l_waveSpeedL,
  l_waveSpeedR);

  REQUIRE(l_waveSpeedL == Approx(-9.731109400863979));
  REQUIRE(l_waveSpeedR == Approx(9.573105166925635));
}

TEST_CASE( "Test the derivation of the fwave net-updates.", "[fwaveUpdates]" ) {

  float l_netUpdatesL[2] = { -5, 3 };
  float l_netUpdatesR[2] = {  4, 7 };
  float l_speed = 0;

  /*
   * Test case (trivial steady state):
   *
   *     left | right
   *   h:  10 | 10
   *  hu:   0 |  0
   */
  tsunami_lab::solvers::fwave::netUpdates( 10,
                                         10,
                                         0,
                                         0,
                                         l_netUpdatesL,
                                         l_netUpdatesR, l_speed);

  REQUIRE( l_netUpdatesL[0] == Approx(0) );
  REQUIRE( l_netUpdatesL[1] == Approx(0) );

  REQUIRE( l_netUpdatesR[0] == Approx(0) );
  REQUIRE( l_netUpdatesR[1] == Approx(0) );
}
