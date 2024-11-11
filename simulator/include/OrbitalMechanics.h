#ifndef ORBITALMECHANICS_H
#define ORBITALMECHANICS_H

#include <array>
#include <vector>
#include <chrono>
#include <Eigen/Dense>

namespace OrbitalMechanics
{

  // Function to convert Keplerian elements to Cartesian coordinates in geocentric-equatorial reference
  std::pair<Eigen::VectorXd, Eigen::VectorXd> keplerian2ijk(
      double sma, // Semi-major axis (Km)
      double ecc, // Eccentricity
      double inc, // Inclination (rad)
      double w,   // Argument of Perigee (rad)
      double nu,  // True Anomaly (rad)
      double raan // Right Ascension of Ascending Node (rad)
  );

  // Calculate eccentric anomalies for every timestamp (in seconds) based on a single mean anomaly, mean motion, eccentricity, and epoch
  Eigen::VectorXd eccentricAnomaly(const std::vector<std::chrono::system_clock::time_point> &timestamps, double M, double mean_motion, double ecc, std::chrono::system_clock::time_point epoch);

  Eigen::VectorXd trueAnomaly(const Eigen::VectorXd &E, double e);
}

#endif
