#include <catch2/catch.hpp>
#define private public
#include "fwave.h"
#undef public

TEST_CASE("Test Eigenvalues lambda", "[lambdas]") {
  /*
   * Test case:
   *  h: 10 | 9
   *  u: -3 | 3
   *
   * roe height: 9.5
   * roe velocity: (sqrt(10) * -3 + 3 * 3) / ( sqrt(10) + sqrt(9) )
   *               = -0.0790021169691720
   * roe speeds: s1 = -0.079002116969172024 - sqrt(9.80665 * 9.5) = -9.7311093998375095
   *             s2 = -0.079002116969172024 + sqrt(9.80665 * 9.5) =  9.5731051658991654
   */

  float l_waveSpeedL = 0;
  float l_waveSpeedR = 0;
  tsunami_lab::solvers::fwave::waveSpeeds(10, 9, -3, 3, l_waveSpeedL,
  l_waveSpeedR);

  REQUIRE( l_waveSpeedL == Approx( -9.7311093998375095 ) );
  REQUIRE( l_waveSpeedR == Approx(  9.5731051658991654 ) );
}
/*
TEST_CASE( "Test the derivation of the fwave wave speeds.", "[fwaveStrengths]" ) {

  float l_strengthL = 0;
  float l_strengthR = 0;

  tsunami_lab::solvers::fwave::waveStrengths( 10,
                                            9,
                                            30,
                                            -27,
                                            -9.7311093998375095,
                                            9.5731051658991654,
                                            l_strengthL,
                                            l_strengthR );

  REQUIRE( l_strengthL == Approx( 41.950951524007173) );
  REQUIRE( l_strengthR == Approx( 15.049048475992827) );
}
*/
TEST_CASE( "Test the derivation of the fwave net-updates.", "[fwaveUpdates]" ) {

  float l_netUpdatesL[2] = { -5, 3 };
  float l_netUpdatesR[2] = {  4, 7 };
  float l_speed = 0;

  /*
  * Test case:
  *  h: 10 | 9
  *  u: -3 | 3

  * d_f1: hur - hul = 27+30= 57
  * d_f2: ((27)²+1/2*9.80665*9²)-((-30)²+1/2*9.80665*10²) = (27)²-(-30)²+1/2*9.80665*(9²-10²) =-264.16317500000000

  * x= 1/(l2-l1) = 1/(9.5731051658991654+9.7311093998375095) = 0.051802159398648326
  * a1 = x * (l2*d_f1-d_f2) = 0.051802159398648326*(9.5731051658991654*(57)+264.16317500000000) = 41.950951524007173
  * a2 = x* (-l1*d_f1+d_f2) = 0.051802159398648326*(9.7311093998375095*(57)-264.16317500000000) = 15.049048475992827

  * z11 = a1
  * z21 = a2
  * z12 = a1 *l1 = 41.950951524007173* (-9.7311093998375095) = -408.22929870739390
  * z22 = a2*l2  = 15.049048475992827* 9.5731051658991654 = 144.06612370739389
  */

  tsunami_lab::solvers::fwave::netUpdates( 10,
                                         9,
                                         -30,
                                         27,
                                         l_netUpdatesL,
                                         l_netUpdatesR,  l_speed );

  REQUIRE( l_netUpdatesL[0] == Approx( 33.5590017014261447899292 ) );
  REQUIRE( l_netUpdatesL[1] == Approx( -326.56631690591093200508 ) );
  REQUIRE( l_netUpdatesR[0] == Approx( 23.4409982985738561366777 ) );
  REQUIRE( l_netUpdatesR[1] == Approx( 224.403141905910928927533 ) );

  /*
  * Test case:
  *  h: 10 | 8
  *  u: 0 | 0
  */

  tsunami_lab::solvers::fwave::netUpdates( 10,
                                         8,
                                         0,
                                         0,
                                         l_netUpdatesL,
                                         l_netUpdatesR, l_speed );

  REQUIRE( l_netUpdatesL[0] ==  Approx(9.394671362) );
  REQUIRE( l_netUpdatesL[1] == -Approx(88.25985)    );

  REQUIRE( l_netUpdatesR[0] == -Approx(9.394671362) );
  REQUIRE( l_netUpdatesR[1] == -Approx(88.25985)    );
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
