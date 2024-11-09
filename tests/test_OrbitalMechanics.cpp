#include <iostream>
#include <array>
#include <cmath>
#include "OrbitalMechanics.h"

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

int main()
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
