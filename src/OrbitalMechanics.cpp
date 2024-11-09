#include "OrbitalMechanics.h"
#include <cmath>

namespace OrbitalMechanics
{
  const double mu_earth = 3.986e5; // Earth's gravitational constant (Km^3/s^2)

  std::pair<std::array<double, 3>, std::array<double, 3>> keplerian2ijk(
      double sma, // Semi-major axis (Km)
      double ecc, // Eccentricity
      double inc, // Inclination (rad)
      double w,   // Argument of Perigee (rad)
      double nu,  // True Anomaly (rad)
      double raan // Right Ascension of Ascending Node (rad)
  )
  {
    // Intermediate calculations
    double a = sma;
    double p = a * (1 - std::pow(ecc, 2)); // Semi-latus rectum
    double r_0 = p / (1 + ecc * cos(nu));  // Distance in the perifocal frame

    // Position components in perifocal reference system (Oxyz)
    double x = r_0 * cos(nu);
    double y = r_0 * sin(nu);

    // Velocity components in perifocal reference system (Oxyz)
    double Vx_ = -std::sqrt(mu_earth / p) * sin(nu);
    double Vy_ = std::sqrt(mu_earth / p) * (ecc + cos(nu));

    // Position components in the geocentric-equatorial reference system (OXYZ)
    double X = (cos(raan) * cos(w) - sin(raan) * sin(w) * cos(inc)) * x +
               (-cos(raan) * sin(w) - sin(raan) * cos(w) * cos(inc)) * y;
    double Y = (sin(raan) * cos(w) + cos(raan) * sin(w) * cos(inc)) * x +
               (-sin(raan) * sin(w) + cos(raan) * cos(w) * cos(inc)) * y;
    double Z = (sin(w) * sin(inc)) * x + (cos(w) * sin(inc)) * y;

    // Velocity components in the geocentric-equatorial reference system (OXYZ)
    double Vx = (cos(raan) * cos(w) - sin(raan) * sin(w) * cos(inc)) * Vx_ +
                (-cos(raan) * sin(w) - sin(raan) * cos(w) * cos(inc)) * Vy_;
    double Vy = (sin(raan) * cos(w) + cos(raan) * sin(w) * cos(inc)) * Vx_ +
                (-sin(raan) * sin(w) + cos(raan) * cos(w) * cos(inc)) * Vy_;
    double Vz = (sin(w) * sin(inc)) * Vx_ + (cos(w) * sin(inc)) * Vy_;

    // Return the results as a pair of arrays
    std::array<double, 3> r = {X, Y, Z};
    std::array<double, 3> v = {Vx, Vy, Vz};

    return {r, v};
  }

} // namespace OrbitalMechanics
