/**
 * @author Julius Isken, Max Engel
 *
 * @section LICENSE
 * Copyright 2020, Julius Isken, Max Engel
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
 * Two-dimensional wave propagation patch.
 **/
#include <catch2/catch.hpp>

#include "WavePropagation2d.h"

TEST_CASE("Test the 2d wave propagation solver in a steady State.",
          "[WaveProp2dSteady]") {
  /*
   * Test case:
   *
   *   Single dam break problem between cell 49 and 50.
   *     left | right
   *       10 | 8
   *        0 | 0
   *
   *   Elsewhere steady state.
   *
   * The net-updates at the respective edge are given as
   * (see derivation in Roe solver):
   *    left          | right
   *      9.394671362 | -9.394671362
   *    -88.25985     | -88.25985
   */

  // construct solver and setup a dambreak problem
  std::size_t l_xRows = 8;
  std::size_t l_yColums = 3;
  tsunami_lab::patches::WavePropagation2d m_waveProp(l_xRows, l_yColums);
  for (std::size_t l_ceY = 0; l_ceY < l_yColums; l_ceY++) {
    for (std::size_t l_ceX = 0; l_ceX < l_xRows; l_ceX++) {
      m_waveProp.setHeight(l_ceX, l_ceY, 10);
      m_waveProp.setMomentumX(l_ceX, l_ceY, 0);
      m_waveProp.setMomentumY(l_ceX, l_ceY, 0);
    }
  }
  

  // perform a time step
  m_waveProp.timeStep(0.1, 1);

  for (std::size_t l_ceY = 0; l_ceY < l_yColums; l_ceY++) {
    // steady state
    for (std::size_t l_ce = 0; l_ce < l_xRows; l_ce++) {
      REQUIRE(m_waveProp.getHeight()[l_ce + l_xRows * l_ceY] == Approx(10));
      REQUIRE(m_waveProp.getMomentumX()[l_ce + l_xRows * l_ceY] == Approx(0));
      REQUIRE(m_waveProp.getMomentumY()[l_ce + l_xRows * l_ceY] == Approx(0));
    }
  }
  /*
    for(std::size_t l_ceY = 0; l_ceY<10; l_ceY++){
        // steady state
      for (std::size_t l_ce = 0; l_ce < 49; l_ce++) {
        REQUIRE(m_waveProp.getHeight()[l_ce + 101 * l_ceY] == Approx(10));
        REQUIRE(m_waveProp.getMomentumX()[l_ce + 101 * l_ceY] == Approx(0));
      }

      // dam-break
      REQUIRE(m_waveProp.getHeight()[49 + 101 * l_ceY] == Approx(10 - 0.1
    * 9.394671362)); REQUIRE(m_waveProp.getMomentumX()[49 + 101 * l_ceY] ==
    Approx(0 + 0.1 * 88.25985));

      REQUIRE(m_waveProp.getHeight()[50 + 101 * l_ceY] == Approx(8 + 0.1
    * 9.394671362)); REQUIRE(m_waveProp.getMomentumX()[50 + 101 * l_ceY] ==
    Approx(0 + 0.1 * 88.25985));

      // steady state
      for (std::size_t l_ce = 51; l_ce < 100; l_ce++) {
        REQUIRE(m_waveProp.getHeight()[l_ce + 101 * l_ceY] == Approx(8));
        REQUIRE(m_waveProp.getMomentumX()[l_ce + 101 * l_ceY] == Approx(0));
      }
    }
    */
}
