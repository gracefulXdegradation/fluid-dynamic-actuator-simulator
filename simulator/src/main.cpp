#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "TLE.h"
#include "config.hpp"
#include "dateTime.hpp"
#include "OrbitalMechanics.h"
#include "DB.hpp"

using namespace Eigen;
using namespace std;

Eigen::VectorXd lla2ecef(const Eigen::VectorXd &lla)
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

    Eigen::VectorXd ecef(3);
    ecef << x, y, z;
    return ecef;
}

MatrixXd cross(const MatrixXd &a, const MatrixXd &b)
{
    if (a.rows() != 3 || b.rows() != 3 || a.cols() != b.cols())
    {
        std::cerr << "Error: Matrices must be 3xN in size and have the same number of columns.\n";
        exit(1); // Exit if the matrices have incompatible sizes
    }

    // Initialize the result matrix with the same number of rows and columns
    MatrixXd x(3, a.cols());

    // Loop over each column and compute the cross product
    for (int i = 0; i < a.cols(); ++i)
    {
        // Extract the 3D vectors from the matrices (column i)
        Vector3d posVec = a.col(i);
        Vector3d velVec = b.col(i);

        // Compute the cross product and store it in the result matrix
        x.col(i) = posVec.cross(velVec);
    }

    return x;
}

std::vector<Matrix3d> reshapeAndPagetranspose(const MatrixXd &x, const MatrixXd &y, const MatrixXd &z)
{
    int cols = x.cols();
    int rows = x.rows();

    std::vector<MatrixXd> reshaped_x, reshaped_y, reshaped_z;

    for (int i = 0; i < cols; ++i)
    {
        MatrixXd slice_x(rows, 1);
        MatrixXd slice_y(rows, 1);
        MatrixXd slice_z(rows, 1);

        for (int j = 0; j < rows; ++j)
        {
            slice_x(j, 0) = x(j, i);
            slice_y(j, 0) = y(j, i);
            slice_z(j, 0) = z(j, i);
        }

        reshaped_x.push_back(slice_x);
        reshaped_y.push_back(slice_y);
        reshaped_z.push_back(slice_z);
    }

    std::vector<Matrix3d> attitude_matrix_pages;

    for (int i = 0; i < cols; ++i)
    {
        // Combine slices along the third dimension (column)
        MatrixXd attitude_matrix(rows, rows);
        attitude_matrix.col(0) = reshaped_x[i];
        attitude_matrix.col(1) = reshaped_y[i];
        attitude_matrix.col(2) = reshaped_z[i];

        attitude_matrix_pages.push_back(attitude_matrix.transpose()); // Page transpose
    }

    return attitude_matrix_pages;
}

std::pair<std::vector<Quaterniond>, MatrixXd> nadir_frame(const MatrixXd &position, const MatrixXd &velocity)
{
    int count = position.cols();

    MatrixXd angular_momentum = cross(position, velocity);

    // z = -position ./ vecnorm(position);
    MatrixXd z = -(position.array().rowwise() / position.colwise().norm().array());
    // y = -angular_momentum ./ vecnorm(angular_momentum);
    MatrixXd y = -(angular_momentum.array().rowwise() / angular_momentum.colwise().norm().array());
    MatrixXd x = cross(y, z);

    std::vector<Matrix3d> attitude_matrix_pages = reshapeAndPagetranspose(x, y, z);

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

void target_pointing_frame(const MatrixXd &position, const MatrixXd &velocity, const Eigen::VectorXd &gs_r, const std::vector<std::chrono::_V2::system_clock::time_point> &date_times)
{
    auto ecef = lla2ecef(gs_r); // in meters

    std::cout << "ECEF Coordinates:\n";
    std::cout << ecef << "\n";
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
