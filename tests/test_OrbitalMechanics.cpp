#include <iostream>
#include <array>
#include <cmath>
#include <chrono>
#include "OrbitalMechanics.h"
#include "dateTime.hpp"

bool areApproxEqual(const std::array<double, 3> &arr1, const std::array<double, 3> &arr2, double tolerance = 1e-6)
{
  for (size_t i = 0; i < arr1.size(); ++i)
  {
    if (std::fabs(arr1[i] - arr2[i]) > tolerance)
    {
      return false;
    }
  }
  return true;
}

bool test_keplerian2ijk()
{
  // Test input values
  double sma = 6852.6;  // Semi-major axis in Km
  double ecc = 0.0013;  // Eccentricity
  double inc = 1.7047;  // Inclination in radians
  double w = 4.7190;    // Argument of Perigee in radians
  double nu = 1.7497;   // True Anomaly in radians
  double raan = 2.7630; // Right Ascension of Ascending Node in radians

  // Expected output (sample values; adjust as necessary)
  std::array<double, 3> expected_position = {-6197.14, 2646.76, 1252.95}; // Example position vector in Km
  std::array<double, 3> expected_velocity = {1.66783, 0.413577, 7.42888}; // Example velocity vector in Km/s

  auto [computed_position, computed_velocity] = OrbitalMechanics::keplerian2ijk(sma, ecc, inc, w, nu, raan);

  // Check if the computed values match the expected values within a tolerance
  if (areApproxEqual(computed_position, expected_position, 1e-2) && areApproxEqual(computed_velocity, expected_velocity, 1e-5))
  {
    std::cout << "Test passed: Computed values match expected values.\n";
  }
  else
  {
    std::cout << "Test failed: Computed values do not match expected values.\n";
    std::cout << "Computed Position: [" << computed_position[0] << ", " << computed_position[1] << ", " << computed_position[2] << "]\n";
    std::cout << "Expected Position: [" << expected_position[0] << ", " << expected_position[1] << ", " << expected_position[2] << "]\n";
    std::cout << "Computed Velocity: [" << computed_velocity[0] << ", " << computed_velocity[1] << ", " << computed_velocity[2] << "]\n";
    std::cout << "Expected Velocity: [" << expected_velocity[0] << ", " << expected_velocity[1] << ", " << expected_velocity[2] << "]\n";
    return 1;
  }

  return 0;
}

bool test_eccentricAnomaly()
{
  double tolerance = 1e-3;

  // Define test parameters
  double initial_mean_anomaly = 1.5638411897256992; // Example mean anomaly (radians)
  double eccentricity = 0.0012642;                  // Example eccentricity
  double mean_motion = 0.0011129785881891408;       // Mean motion (radians per second)
  auto epoch = DateTime::parseDateTime("2023-06-26 08:46:57");

  auto startTime = DateTime::parseDateTime("2023-07-04 14:25:00");
  int timeStep = 500;

  auto dur = std::chrono::duration_cast<std::chrono::seconds>(startTime - epoch).count();

  std::vector<std::chrono::system_clock::time_point> timestamps = {
      startTime,
      startTime + std::chrono::milliseconds(timeStep),
      startTime + std::chrono::milliseconds(timeStep * 2),
      startTime + std::chrono::milliseconds(timeStep * 3),
      startTime + std::chrono::milliseconds(timeStep * 4)};

  std::vector<double> expected = {
      1.748486401596818,
      1.749042766375569,
      1.749599130769283,
      1.750155494778001,
      1.750711858401648};

  // Calculate the eccentric anomalies
  auto eccentricAnomalies = OrbitalMechanics::eccentricAnomaly(timestamps, initial_mean_anomaly, mean_motion, eccentricity, epoch);

  bool isFail = false;
  // Check results
  for (size_t i = 0; i < eccentricAnomalies.size(); ++i)
  {
    if (std::fabs(eccentricAnomalies[i] - expected[i]) > tolerance)
    {
      isFail = true;
    }
  }

  if (isFail)
  {
    std::cout << "Test failed: Computed values do not match expected values.\n";
    for (size_t i = 0; i < eccentricAnomalies.size(); ++i)
    {
      std::cout << "Eccentric Anomaly at t[" << i << "]: " << eccentricAnomalies[i] << " rad. Expected: " << expected[i] << "\n";
    }
  }
  else
  {
    std::cout << "Test passed: Computed values match expected values.\n";
  }

  return isFail;
}

int main()
{
  bool result_keplerian2ijk = test_keplerian2ijk();
  bool result_eccentricAnomaly = test_eccentricAnomaly();
  return result_keplerian2ijk || result_eccentricAnomaly;
}
