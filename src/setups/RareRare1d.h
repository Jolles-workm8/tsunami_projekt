#ifndef TSUNAMI_LAB_SETUPS_RARE_RARE_1D_H
#define TSUNAMI_LAB_SETUPS_RARE_RARE_1D_H

#include "Setup.h"

namespace tsunami_lab {
namespace setups {
class RareRare1d;
}
} // namespace tsunami_lab

/**
 * 1d Rarerare problem
 **/

class tsunami_lab::setups::RareRare1d : public Setup {
private:
  //! impuls of the wave
  t_real m_impuls = 0;

  //! location of the impuls
  t_real m_location = 0;

  //! water height in x
  t_real m_height = 0;

public:
  /**
   *Constructor
   * @param i_impuls impuls of the waves namely hu
   * @param i_location position where the waves start
   * @param i_height water height
   **/
  RareRare1d(t_real i_impuls, t_real i_location, t_real i_height);

  /**
   * Gets the water height at a given point.
   *
   * @param i_x x-coordinate of the queried point.
   * @return height at the given point.
   **/
  t_real getHeight(t_real, t_real) const;

  /**
   * Gets the momentum in x-direction.
   *
   * @return momentum in x-direction.
   **/
  t_real getMomentumX(t_real i_x, t_real) const;

  /**
   * Gets the momentum in y-direction.
   *
   * @return momentum in y-direction.
   **/
  t_real getMomentumY(t_real, t_real) const;
};

#endif
