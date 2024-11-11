#include <iostream>
#include <string>
#include <Eigen/Dense>
#include "TLE.h"
#include "config.hpp"
#include "dateTime.hpp"
#include "OrbitalMechanics.h"
#include "DB.hpp"

using namespace Eigen;
using namespace std;

Eigen::MatrixXd cross(const Eigen::MatrixXd &a, const Eigen::MatrixXd &b)
{
    if (a.rows() != 3 || b.rows() != 3 || a.cols() != b.cols())
    {
        std::cerr << "Error: Matrices must be 3xN in size and have the same number of columns.\n";
        exit(1); // Exit if the matrices have incompatible sizes
    }

    // Initialize the result matrix with the same number of rows and columns
    Eigen::MatrixXd x(3, a.cols());

    // Loop over each column and compute the cross product
    for (int i = 0; i < a.cols(); ++i)
    {
        // Extract the 3D vectors from the matrices (column i)
        Eigen::Vector3d posVec = a.col(i);
        Eigen::Vector3d velVec = b.col(i);

        // Compute the cross product and store it in the result matrix
        x.col(i) = posVec.cross(velVec);
    }

    return x;
}

void nadir_frame(const Eigen::MatrixXd &position, const Eigen::MatrixXd &velocity)
{
    Eigen::MatrixXd angular_momentum = cross(position, velocity);

    // z = -position ./ vecnorm(position);
    Eigen::MatrixXd z = -(position.array().rowwise() / position.colwise().norm().array());
    // y = -angular_momentum ./ vecnorm(angular_momentum);
    Eigen::MatrixXd y = -(angular_momentum.array().rowwise() / angular_momentum.colwise().norm().array());
    Eigen::MatrixXd x = cross(y, z);

    std::cout << x.col(0) << std::endl;
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

        Eigen::MatrixXd m_i_r(3, date_times.size());
        m_i_r.setZero();
        Eigen::MatrixXd m_i_v(3, date_times.size());
        m_i_v.setZero();

        for (int i = 0; i < trueAnomalies.size(); i++)
        {
            auto &ta = trueAnomalies[i];
            auto [i_r, i_v] = OrbitalMechanics::keplerian2ijk(tle.getSemiMajorAxis(), tle.getEccentricity(), tle.getInclination(), tle.getArgumentOfPerigee(), ta, tle.getRightAscension());
            m_i_r.col(i) = i_r;
            m_i_v.col(i) = i_v;
        }

        nadir_frame(m_i_r, m_i_v);

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
