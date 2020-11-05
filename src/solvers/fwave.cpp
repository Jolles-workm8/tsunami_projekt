#include "fwave.h"
#include <cmath>

void tsunami_lab::solvers::fwave::lambdaRoe(t_real ul, t_real ur, t_real hr,
                                            t_real hl, &t_real lambda) {

  // precompute square roots
  t_real hl_sqrt = std::sqrt(hl);
  t_real hr_sqrt = std::sqrt(hr);

  // precompute h_roe and u_roe

  t_real h_roe_sqrt = std::sqrt((hl + hr) / 2);
  t_real u_roe = (ul * hl_sqrt + ur * hr_sqrt) / (hr_sqrt + hl_sqrt);

  // compute the lambdas

  lambda[0] = u_roe - (m_gSqrt * h_roe_sqrt);
  lambda[1] = u_roe + (m_gSqrt * h_roe_sqrt);
}

void tsunami_lab::solvers::fwave::matrices(t_real ul, t_real ur, t_real hr,
                                           t_real hl) {

  // compute lambdas
  t_real lambda[2];
}
