#ifndef ORBITALMECHANICS_H
#define ORBITALMECHANICS_H

#include <array>

namespace OrbitalMechanics
{

  // Function to convert Keplerian elements to Cartesian coordinates in geocentric-equatorial reference
  std::pair<std::array<double, 3>, std::array<double, 3>> keplerian2ijk(
      double sma, // Semi-major axis (Km)
      double ecc, // Eccentricity
      double inc, // Inclination (rad)
      double w,   // Argument of Perigee (rad)
      double nu,  // True Anomaly (rad)
      double raan // Right Ascension of Ascending Node (rad)
  );

} // namespace OrbitalMechanics

#endif // ORBITALMECHANICS_H
