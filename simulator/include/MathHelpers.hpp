#ifndef MATHHELPERS_HPP
#define MATHHELPERS_HPP

#include <vector>
#include <cmath>

namespace MathHelpers
{
  inline double wrapTo2Pi(double angle)
  {
    return fmod(angle + 2 * M_PI, 2 * M_PI);
  }

  inline double deg2rad(double deg)
  {
    return deg * M_PI / 180.0;
  }

  inline double rad2deg(double rad)
  {
    return rad * 180.0 / M_PI;
  }

}

#endif
