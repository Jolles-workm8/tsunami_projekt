#include "fwave.h"
#include <cmath>

void tsunami_lab::solvers::fwave::zsolver(t_real i_hL, t_real i_hR,
                                          t_real i_huL, t_real i_huR, t_real &z,
                                          t_real &lambda) {

  // prepare the variables
  t_real i_uL = i_huL / i_hL;
  t_real i_uR = i_huR / i_hR;

  // precompute square roots
  t_real hl_sqrt = std::sqrt(i_hL);
  t_real hr_sqrt = std::sqrt(i_hR);

  // precompute h_roe and u_roe
  t_real h_roe_sqrt = std::sqrt((i_hL + i_hR) / 2);
  t_real u_roe = (i_uL * hl_sqrt + i_uR * hr_sqrt) / (hr_sqrt + hl_sqrt);

  // compute the lambdas using array to align the memory for parallelization
  // later
  lambda[0] = u_roe - (m_g Sqrt * h_roe_sqrt);
  lambda[1] = u_roe + (m_gSqrt * h_roe_sqrt);

  // precompute the factors that can be precomputed, hard square is faster then
  // casting power in C++
  t_real lambda_factor = 1 / (lambda[1] - lambda[0]);
  t_real df1 = i_huR - i_huL;
  t_real df2 = i_huR * i_huR - i_huL * i_huL + g / 2 * (i_hR - i_hL);

  // compute the z's, again array to align memory, for formula see documentation
  z[0] = lambda_factor * ((lambda[1] * df1) - lambda[0] * lambda[0] * df1 +
                          (lambda[0] - 1) * df2);
  z[1] = lambda_factor * ((lambda[1] * df1) - lambda[1] * lambda[0] * df1 +
                          (lambda[1] - 1) * df2);
}

void tsunami_lab::solvers::fwave::netUpdates(t_real i_hL, t_real i_hR,
                                             t_real i_huL, t_real i_huR,
                                             t_real o_netUpdateL[2],
                                             t_real o_netUpdateR[2]) {

  // set memory for lambda and z.

  t_real lambda[2];
  t_real z[2];
  // compute the z's

  zsolver(i_hL, i_hR, i_huL, i_huR, z, lambda);

  for (int i = 0, i < 2, i++) {
    if (lambda[i] < 0) {
      a_dQ_minus += z[i];
    } else {
      a_dq_plus += z[i];
    }
  }
}
