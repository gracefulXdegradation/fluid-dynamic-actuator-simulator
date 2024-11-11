#include <cmath>
#include "OrbitalMechanics.h"
#include "MathHelpers.hpp"
#include <Eigen/Dense>

namespace OrbitalMechanics
{
  const double mu_earth = 3.986e5; // Earth's gravitational constant (Km^3/s^2)

  std::pair<Eigen::VectorXd, Eigen::VectorXd> keplerian2ijk(
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
    // std::array<double, 3> r = {X, Y, Z};
    // std::array<double, 3> v = {Vx, Vy, Vz};

    Eigen::VectorXd r(3);
    r << X, Y, Z;
    Eigen::VectorXd v(3);
    v << Vx, Vy, Vz;

    return {r, v};
  }

  Eigen::VectorXd eccentricAnomaly(const std::vector<std::chrono::system_clock::time_point> &timestamps, double M, double mean_motion, double ecc, std::chrono::system_clock::time_point epoch)
  {
    const double tolerance = 1e-9;
    Eigen::VectorXd eccentricAnomalies(timestamps.size());
    eccentricAnomalies.setZero();

    for (int i = 0; i < timestamps.size(); i++)
    {
      double elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(timestamps[i] - epoch).count() / 1000.0;

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

      eccentricAnomalies[i] = MathHelpers::wrapTo2Pi(Ens);
    }

    return eccentricAnomalies;
  }

  Eigen::VectorXd trueAnomaly(const Eigen::VectorXd &E, double e)
  {
    Eigen::VectorXd nu(E.size());
    nu.setZero();

    for (int i = 0; i < E.size(); i++)
    {
      auto eccentric_anomaly = E[i];
      double sin_part = std::sqrt(1 - e * e) * std::sin(eccentric_anomaly);
      double cos_part = std::cos(eccentric_anomaly) - e;
      double true_anom = std::atan2(sin_part, cos_part);

      nu[i] = MathHelpers::wrapTo2Pi(true_anom);
    }

    return nu;
  }
}
