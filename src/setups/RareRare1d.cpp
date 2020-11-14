#include "RareRare1d.h"

tsunami_lab::setups::RareRare1d::RareRare1d(t_real i_impuls, t_real i_location,
                                            t_real i_height) {
  m_impuls = i_impuls;
  m_location = i_location;
  m_height = i_height;
}

tsunami_lab::t_real tsunami_lab::setups::RareRare1d::getHeight(t_real,
                                                               t_real) const {
  return m_height;
}

tsunami_lab::t_real
tsunami_lab::setups::RareRare1d::getMomentumX(t_real i_x, t_real) const {
  if (i_x < m_location) {
    return -m_impuls;
  } else {
    return m_impuls;
  }
}

tsunami_lab::t_real
tsunami_lab::setups::RareRare1d::getMomentumY(t_real, t_real) const {
  return 0;
}
