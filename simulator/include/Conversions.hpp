#ifndef CONVERSIONS_HPP
#define CONVERSIONS_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

namespace Conversions
{
  // Constants with high precision
  const double ARC_SECONDS_TO_RADIANS = M_PI / 648000.0;
  const double EARTH_ROTATION_DERIVATIVE = M_PI * 1.00273781191135448 / 43200.0;
  const Matrix3d DERIVATIVE_MATRIX = (Matrix3d() << 0.0, -EARTH_ROTATION_DERIVATIVE, 0.0,
                                      EARTH_ROTATION_DERIVATIVE, 0.0, 0.0,
                                      0.0, 0.0, 0.0)
                                         .finished();

  // Helper: Julian Date calculation
  inline double utc_time_to_julian_date(std::chrono::system_clock::time_point utc_time)
  {
    auto time_t = std::chrono::system_clock::to_time_t(utc_time);
    auto tm = *std::gmtime(&time_t);

    int year = tm.tm_year + 1900;
    int month = tm.tm_mon + 1;
    int day = tm.tm_mday;
    int hour = tm.tm_hour;
    int minute = tm.tm_min;
    int second = tm.tm_sec;

    double julian_date = 367 * year - 7 * (year + (month + 9) / 12) / 4 + 275 * month / 9 + day + 1721013.5;
    julian_date += (hour + minute / 60.0 + second / 3600.0) / 24.0;
    return julian_date;
  }

  // Helper: Polynomial evaluation
  inline double evaluate_polynomial(const std::vector<double> &coeffs, double x)
  {
    double result = 0.0;
    for (int i = coeffs.size() - 1; i >= 0; --i)
    {
      result = result * x + coeffs[i];
    }
    return result;
  }

  // Precession functions
  inline double precession_x(double julian_century)
  {
    return evaluate_polynomial({-0.016617, 2004.191898, -0.4297829, -0.19861834}, julian_century);
  }

  inline double precession_y(double julian_century)
  {
    return evaluate_polynomial({-0.006951, -0.025896, -22.4072747, 0.00190059}, julian_century);
  }

  // Nutation/Precession correction matrix
  inline std::pair<double, double> compute_celestial_positions(double julian_century)
  {
    double celestial_x = precession_x(julian_century) * ARC_SECONDS_TO_RADIANS;
    double celestial_y = precession_y(julian_century) * ARC_SECONDS_TO_RADIANS;
    return {celestial_x, celestial_y};
  }

  // Rotation Matrix
  inline Matrix3d rotation_matrix_ecef_to_eci(double julian_century)
  {
    double earth_rotation_angle = 2 * M_PI * (0.7790572732640 + 1.00273781191135448 * 36525.0 * julian_century);
    Matrix3d earth_matrix = Matrix3d::Identity();
    earth_matrix(0, 0) = earth_matrix(1, 1) = std::cos(earth_rotation_angle);
    earth_matrix(1, 0) = std::sin(earth_rotation_angle);
    earth_matrix(0, 1) = -earth_matrix(1, 0);

    auto [gcrs_x, gcrs_y] = compute_celestial_positions(julian_century);
    double a = 0.5 + 0.125 * (gcrs_x * gcrs_x + gcrs_y * gcrs_y);

    Matrix3d pn_matrix;
    pn_matrix << 1 - a * gcrs_x * gcrs_x, -a * gcrs_x * gcrs_y, gcrs_x,
        -a * gcrs_x * gcrs_y, 1 - a * gcrs_y * gcrs_y, gcrs_y,
        -gcrs_x, -gcrs_y, 1 - a * (gcrs_x * gcrs_x + gcrs_y * gcrs_y);

    return pn_matrix * earth_matrix;
  }

  // Main transformation function
  inline std::pair<Vector3d, Vector3d> ecef_to_eci(
      const Vector3d &ecef_point,
      std::chrono::system_clock::time_point utc_time,
      const Vector3d &ecef_velocity = Vector3d::Zero())
  {
    double julian_day = utc_time_to_julian_date(utc_time);
    double julian_century = (julian_day - 2451545.0) / 36525.0;

    Matrix3d rotation_ecef_to_eci = rotation_matrix_ecef_to_eci(julian_century);
    Vector3d eci_point = rotation_ecef_to_eci * ecef_point;

    // Velocity transformation
    Vector3d eci_velocity = rotation_ecef_to_eci * ecef_velocity + (rotation_ecef_to_eci * DERIVATIVE_MATRIX) * ecef_point;

    return {eci_point, eci_velocity};
  }

  // https://fossies.org/linux/gpsd/contrib/lla2ecef.c
  inline VectorXd lla2ecef(const VectorXd &lla)
  {
    const double lat = lla[0];
    const double lon = lla[1];
    const double alt = lla[2];

    // WGS84 ellipsoid constants
    constexpr double a = 6378137.0;         // semi-major axis in meters
    constexpr double f = 1 / 298.257223563; // flattening
    constexpr double e2 = 2 * f - f * f;    // first eccentricity squared

    // Convert latitude and longitude from degrees to radians
    double latRad = lat * M_PI / 180.0;
    double lonRad = lon * M_PI / 180.0;

    // Radius of curvature in the prime vertical
    double N = a / std::sqrt(1 - e2 * std::sin(latRad) * std::sin(latRad));

    double x = (N + alt) * std::cos(latRad) * std::cos(lonRad);
    double y = (N + alt) * std::cos(latRad) * std::sin(lonRad);
    double z = ((1 - e2) * N + alt) * std::sin(latRad);

    VectorXd ecef(3);
    ecef << x, y, z;
    return ecef;
  }

}

#endif
