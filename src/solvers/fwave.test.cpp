#include <catch2/catch.hpp>
#define private public
#include "fwave.h"
#undef public

TEST_CASE("Test Eigenvalues lambda", "[lambdas]"){
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


  float lambda[2];
  lambda[0] = 0;
  lambda[1] = 0;
  tsunami_lab::solvers::fwave::getLambda(10,
                                         9,
                                         -3,
                                         3,
                                         lambda);

  REQUIRE( lambda[0] == Approx( -9.731109400863979 ) );
  REQUIRE( lambda[1] == Approx(  9.573105166925635 ) );

}


TEST_CASE("Test steady states", "[steadyState]"){

}
