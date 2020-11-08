#include "fwave.h"
#include <cmath>

void tsunami_lab::solvers::fwave::zsolver(t_real &ql, t_real &qr, t_real &z,
                                          t_real &lambda) {

  // TODO: alles neu in einer struct einbauen um besser memory zu loeschen
  // prepare the variables
  t_real hl = ql[0];
  t_real hr = qr[0];

  t_real hul = ql[1];
  t_real hur = qr[1];

  t_real ul = ql[1] / ql[0];
  t_real ur = qr[1] / qr[0];

  // precompute square roots
  t_real hl_sqrt = std::sqrt(hl);
  t_real hr_sqrt = std::sqrt(hr);

  // precompute h_roe and u_roe

  t_real h_roe_sqrt = std::sqrt((hl + hr) / 2);
  t_real u_roe = (ul * hl_sqrt + ur * hr_sqrt) / (hr_sqrt + hl_sqrt);

  // compute the lambdas using array to align the memory for parallelization
  // later
  lambda[0] = u_roe - (m_gSqrt * h_roe_sqrt);
  lambda[1] = u_roe + (m_gSqrt * h_roe_sqrt);

  // precompute the factors that can be precomputed, hard square is faster then
  // casting power in C++
  t_real lambda_factor = 1 / (lambda[1] - lambda[0]);
  t_real f1 = hur - hul;
  t_real f2 = hur * hur - hul * hul + g / 2 * (hr - hl);

  // compute the z's, again array to align memory, for formula see documentation
  z[0] = lambda_factor *
         ((lambda[1] * f1) - lambda[0] * lambda[0] * f1 + (lambda[0] - 1) * f2);
  z[1] = lambda_factor *
         ((lambda[1] * f1) - lambda[1] * lambda[0] * f1 + (lambda[1] - 1) * f2);
}

void tsunami_lab::solvers::fwave::netUpdates(t_real &ql, t_real &qr,
                                             t_real net_UpdateL,
                                             t_real net_UpdateR) {

  // set memory for lambda
  // TODO: struct machen für besseres memory managment, delete und new siehe cpp
  // doku

  t_real lambda[2];
  t_real z[2];

  // compute the z's

  zsolver(ql, qr, z, lambda);

  for (int i = 0, i < 2, i++) {
    if (lambda[i] < 0) {
      net_UpdateL += z[i];
    } else {
      net_UpdateR += z[i];
    }
  }
}
