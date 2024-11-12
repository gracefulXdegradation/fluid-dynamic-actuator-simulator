#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "TLE.h"
#include "config.hpp"
#include "dateTime.hpp"
#include "OrbitalMechanics.h"
#include "DB.hpp"
#include "MathHelpers.hpp"
#include "Conversions.hpp"

using namespace Eigen;
using namespace std;

std::pair<std::vector<Quaterniond>, MatrixXd> nadir_frame(const MatrixXd &position, const MatrixXd &velocity)
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

void target_pointing_frame(const MatrixXd &r_inert, const MatrixXd &v_inert, const Vector3d &gs_r, const std::vector<std::chrono::_V2::system_clock::time_point> &date_times)
{
    Vector3d gs_r_ecef = Conversions::lla2ecef(gs_r); // in meters
    Vector3d gs_v_ecef(3);
    gs_v_ecef << 0, 0, 0;
    MatrixXd gs_r_inert(3, date_times.size());
    MatrixXd gs_v_inert(3, date_times.size());
    for (int i = 0; i < date_times.size(); i++)
    {
        auto [gs_r_eci, gs_v_eci] = Conversions::ecef_to_eci(gs_r_ecef, date_times[i], gs_v_ecef);
        gs_r_inert.col(i) = gs_r_eci;
        gs_v_inert.col(i) = gs_v_eci;
    }

    MatrixXd orbit_angular_momentum = MathHelpers::cross(r_inert, v_inert);

    // ground station position and velocity are converted to km and km/s
    MatrixXd delta_v_inert = gs_v_inert / 1000.0 - v_inert;
    MatrixXd distance_inert = gs_r_inert / 1000.0 - r_inert;

    MatrixXd inertial_target_rate(3, date_times.size());

    for (int i = 0; i < date_times.size(); i++)
    {
        Vector3d dist_inert_vec = distance_inert.col(i).head<3>();
        Vector3d delta_v_inert_vec = delta_v_inert.col(i).head<3>();
        inertial_target_rate.col(i) = dist_inert_vec.cross(delta_v_inert_vec) / dist_inert_vec.dot(dist_inert_vec);
    }

    std::cout << "distance_inert:" << std::endl;
    std::cout << distance_inert.col(0) << std::endl;
    std::cout << "delta_v_inert:" << std::endl;
    std::cout << delta_v_inert.col(0) << std::endl;
    std::cout << "inertial_target_rate:" << std::endl;
    std::cout << inertial_target_rate.col(0) << std::endl;
}

int main()
{
    // Create an instance of ConfigParser with the path to your config.json
    Config config(string(BUILD_OUTPUT_PATH) + "/config.json");

    // Generate discrete time points using the function
    auto date_times = DateTime::generateTimePoints(config.getStartDateTime(), config.getEndDateTime(), config.getControlTimeStep());

    // Output the time points
    cout << "Executing simulation" << endl;
    cout << "====================" << endl;
    cout << "Start date: " << DateTime::formatTime(date_times.at(0)) << endl;

    try
    {
        TLE tle = TLE::fromFile(string(BUILD_OUTPUT_PATH) + "/tle.txt");

        auto eccentricAnomalies = OrbitalMechanics::eccentricAnomaly(date_times, tle.getMeanAnomaly(), tle.getMeanMotion(), tle.getEccentricity(), tle.getEpoch());
        auto trueAnomalies = OrbitalMechanics::trueAnomaly(eccentricAnomalies, tle.getEccentricity());

        MatrixXd m_i_r(3, date_times.size());
        m_i_r.setZero();
        MatrixXd m_i_v(3, date_times.size());
        m_i_v.setZero();

        for (int i = 0; i < trueAnomalies.size(); i++)
        {
            auto &ta = trueAnomalies[i];
            auto [i_r, i_v] = OrbitalMechanics::keplerian2ijk(tle.getSemiMajorAxis(), tle.getEccentricity(), tle.getInclination(), tle.getArgumentOfPerigee(), ta, tle.getRightAscension());
            m_i_r.col(i) = i_r;
            m_i_v.col(i) = i_v;
        }

        auto [q_in, n_omega_n] = nadir_frame(m_i_r, m_i_v);
        target_pointing_frame(m_i_r, m_i_v, config.getGroundStationPosition(), date_times);

        // Save to file
        auto ts = DateTime::getCurrentTimestamp();
        DB::writematrix(m_i_r, "./output/" + ts, "i_r.csv");
        DB::writematrix(m_i_v, "./output/" + ts, "i_v.csv");
        DB::serializeTimePointsToCSV(date_times, "./output/" + ts, "t.csv");

        std::cout << "Data saved to output.txt" << std::endl;
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
