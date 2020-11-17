#include "DamBreakNew.h"

tsunami_lab::setups::DamBreakNew::DamBreakNew(t_real i_heightLeft,
                                              t_real i_heightRight,
                                              t_real i_locationDam,
                                              t_real i_momentumRiver) {
  m_heightLeft = i_heightLeft;
  m_heightRight = i_heightRight;
  m_locationDam = i_locationDam;
  m_momentumRiver = i_momentumRiver;
}

tsunami_lab::t_real tsunami_lab::setups::DamBreakNew::getHeight(t_real i_x,
                                                                t_real) const {
  if (i_x < m_locationDam) {
    return m_heightLeft;
  } else {
    return m_heightRight;
  }
}

tsunami_lab::t_real
tsunami_lab::setups::DamBreakNew::getMomentumX(t_real i_x, t_real) const {
  if (i_x < m_locationDam) {
    return 0;
  } else {
    return m_momentumRiver;
  }
}

tsunami_lab::t_real
tsunami_lab::setups::DamBreakNew::getMomentumY(t_real, t_real) const {
  return 0;
}
