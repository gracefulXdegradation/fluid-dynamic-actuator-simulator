#ifndef FRAMES_HPP
#define FRAMES_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "MathHelpers.hpp"
#include "Conversions.hpp"

using namespace Eigen;
using namespace std;

namespace Frames
{

  inline std::pair<std::vector<Quaterniond>, Matrix3Xd> nadir_frame(const Matrix3Xd &position, const Matrix3Xd &velocity)
  {
    int count = position.cols();

    MatrixXd angular_momentum = MathHelpers::cross(position, velocity);

    // z = -position ./ vecnorm(position);
    MatrixXd z = -(position.array().rowwise() / position.colwise().norm().array());
    // y = -angular_momentum ./ vecnorm(angular_momentum);
    MatrixXd y = -(angular_momentum.array().rowwise() / angular_momentum.colwise().norm().array());
    MatrixXd x = MathHelpers::cross(y, z);

    std::vector<Matrix3d> attitude_matrix_pages = MathHelpers::reshapeAndPagetranspose(x, y, z);

    // attitude = smooth_quaternion(quaternion(attitude_matrix, 'rotmat', 'frame'));
    // ignoring smooth_quaternion
    std::vector<Quaterniond> attitude;
    attitude.reserve(count);
    for (int i = 0; i < count; ++i)
    {
      Matrix3d rot_matrix = attitude_matrix_pages[i];

      // transpose to counter a negative scalar part that gets produced otherwise
      Quaterniond q(rot_matrix.transpose());
      attitude.push_back(q);
    }

    // angular_rate = rotateframe(attitude, (angular_momentum ./ vecnorm(position).^2)')';
    MatrixXd momentum_normalized = angular_momentum.array().rowwise() / position.colwise().norm().array().square().matrix().array();
    MatrixXd angular_rate(angular_momentum.rows(), angular_momentum.cols());
    for (int i = 0; i < count; ++i)
    {
      angular_rate.col(i) = attitude[i].conjugate() * momentum_normalized.col(i);
    }

    return {attitude, angular_rate};
  }

  inline std::tuple<std::vector<Quaterniond>, Matrix3Xd, Matrix3Xd, Matrix3Xd> target_pointing_frame(const MatrixXd &r_inert, const MatrixXd &v_inert, const Vector3d &gs_r, const std::vector<std::chrono::_V2::system_clock::time_point> &date_times)
  {
    Vector3d gs_r_ecef = Conversions::lla2ecef(gs_r); // in meters
    Vector3d gs_v_ecef(3);
    gs_v_ecef << 0, 0, 0;
    Matrix3Xd gs_r_inert(3, date_times.size());
    Matrix3Xd gs_v_inert(3, date_times.size());
    for (int i = 0; i < date_times.size(); i++)
    {
      auto [gs_r_eci, gs_v_eci] = Conversions::ecef_to_eci(gs_r_ecef, date_times[i], gs_v_ecef);
      // ground station position and velocity are converted to km and km/s
      gs_r_inert.col(i) = gs_r_eci / 1000.0;
      gs_v_inert.col(i) = gs_v_eci / 1000.0;
    }

    Matrix3Xd orbit_angular_momentum = MathHelpers::cross(r_inert, v_inert);

    Matrix3Xd delta_v_inert = gs_v_inert - v_inert;
    Matrix3Xd distance_inert = gs_r_inert / 1000.0 - r_inert;

    Matrix3Xd inertial_target_rate(3, date_times.size());

    for (int i = 0; i < date_times.size(); i++)
    {
      Vector3d dist_inert_vec = distance_inert.col(i);
      Vector3d delta_v_inert_vec = delta_v_inert.col(i);
      inertial_target_rate.col(i) = dist_inert_vec.cross(delta_v_inert_vec) / dist_inert_vec.dot(dist_inert_vec);
    }

    Matrix3Xd target_angular_rate(3, date_times.size());
    std::vector<Quaterniond> target_attitude;

    for (int i = 0; i < date_times.size(); i++)
    {
      Vector3d vec = distance_inert.col(i);
      Vector3d ang_mom_norm = MathHelpers::normalizeVector(orbit_angular_momentum.col(i));
      Vector3d izt = MathHelpers::normalizeVector(vec);
      Vector3d iztph = izt.dot(ang_mom_norm) * ang_mom_norm;
      Vector3d iztnh = izt - iztph;
      Vector3d ixt = MathHelpers::normalizeVector(iztnh.cross(orbit_angular_momentum.col(i)));
      Vector3d iyt = MathHelpers::normalizeVector(izt.cross(ixt));

      Eigen::Matrix3d mat;
      mat.col(0) = ixt;
      mat.col(1) = iyt;
      mat.col(2) = izt;
      Quaterniond q(mat);
      target_attitude.push_back(q);
      target_angular_rate.col(i) = q.conjugate() * inertial_target_rate.col(i); // rotate frame
    }

    return std::make_tuple(target_attitude, target_angular_rate, gs_r_inert, gs_v_inert);
  }

}

#endif
