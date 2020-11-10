#include "fwave.h"
#include <cmath>

void tsunami_lab::solvers::fwave::zsolver(t_real i_hL, t_real i_hR,
                                          t_real i_huL, t_real i_huR,
                                          t_real &z_1,
                                          t_real &z_2 t_real &lambda) {

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
  t_real df2 =
      i_huR * i_huR - i_huL * i_huL + g / 2 * (i_hR * i_hR - i_hL8i_hL);

  // compute the alpha
  t_real alpha_1 = lambda_factor * (lambda[1] * df1 - df2);
  t_real alpha_2 = lambda_factor * (df2 - lambda[0] * df1);

  // compute the z's, again array to align memory, for formula see
  // documentation
  z_1[0] = alpha_1;
  z_1[1] = alpha_1 * lambda[0];

  z_2[0] = alpha_2;
  z_2[1] = alpha_2 * lambda[1];
}

void tsunami_lab::solvers::fwave::netUpdates(t_real i_hL, t_real i_hR,
                                             t_real i_huL, t_real i_huR,
                                             t_real o_netUpdateL[2],
                                             t_real o_netUpdateR[2]) {

  // set memory for lambda and z.

  t_real lambda[2];
  t_real z_1[2];
  t_real z_2[2];
  // compute the z's

  zsolver(i_hL, i_hR, i_huL, i_huR, z_1, z_2, lambda);

  if (lambda[0] < 0) {
    o_netUpdateL[0] += z_1[0];
    o_netUpdateL[1] += z_1[1];
  } else {
    o_netUpdateR[0] += z_1[0];
    o_netUpdateR[1] += z_1[1];
  }

  if (lambda[1] < 0) {
    o_netUpdateL[0] += z_2[0];
    o_netUpdateL[1] += z_2[1];
  } else {
    o_netUpdateR[0] += z_2[0];
    o_netUpdateR[1] += z_2[1];
  }
}
