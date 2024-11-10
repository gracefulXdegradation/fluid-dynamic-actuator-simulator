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

        std::vector<std::array<double, 3>> radius_vectors;
        std::vector<std::array<double, 3>> velocity_vectors;
        radius_vectors.reserve(date_times.size());
        velocity_vectors.reserve(date_times.size());

        for (const auto &trueAnomaly : trueAnomalies)
        {
            auto [i_r, i_v] = OrbitalMechanics::keplerian2ijk(tle.getSemiMajorAxis(), tle.getEccentricity(), tle.getInclination(), tle.getArgumentOfPerigee(), trueAnomaly, tle.getRightAscension());
            radius_vectors.push_back(i_r);
            velocity_vectors.push_back(i_v);
        }

        // Save to file
        DB::writematrix(radius_vectors, "./output", "i_r.csv");
        DB::writematrix(velocity_vectors, "./output", "i_v.csv");

        std::cout << "Data saved to output.txt" << std::endl;

        // MatrixXd A(2, 3);
        // A << 1, 2, 3,
        //     4, 5, 6;

        // MatrixXd B(3, 2);
        // B << 7, 8,
        //     9, 10,
        //     11, 12;

        // MatrixXd result = A * B;
        // MatrixXd transposedResult = result.transpose();

        // cout << "Result of (A * B)^T:\n"
        //      << transposedResult << endl;
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
