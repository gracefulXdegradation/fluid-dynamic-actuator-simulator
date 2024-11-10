#include <cmath>
#include "OrbitalMechanics.h"
#include "MathHelpers.hpp"

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

  std::vector<double> eccentricAnomaly(const std::vector<std::chrono::system_clock::time_point> &timestamps, double M, double mean_motion, double ecc, std::chrono::system_clock::time_point epoch)
  {
    const double tolerance = 1e-9;
    std::vector<double> eccentricAnomalies;
    eccentricAnomalies.reserve(timestamps.size());

    for (const auto &timestamp : timestamps)
    {
      double elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp - epoch).count() / 1000.0;

      // Calculate the mean anomaly for the current timestamp
      double current_mean_anomaly = MathHelpers::wrapTo2Pi(M + elapsed_seconds * mean_motion);

      // Solve for the eccentric anomaly using Newton's method
      double En = current_mean_anomaly;
      double Ens = En - (En - ecc * std::sin(En) - current_mean_anomaly) / (1 - ecc * std::cos(En));

      while (std::abs(Ens - En) > tolerance)
      {
        En = Ens;
        Ens = En - (En - ecc * std::sin(En) - current_mean_anomaly) / (1 - ecc * std::cos(En));
      }

      eccentricAnomalies.push_back(MathHelpers::wrapTo2Pi(Ens));
    }

    return eccentricAnomalies;
  }

}
