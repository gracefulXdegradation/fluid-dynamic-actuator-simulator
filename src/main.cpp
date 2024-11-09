#include <iostream>
#include <string>
#include <Eigen/Dense>
#include "TLE.h"
#include "config.hpp"
#include "dateTime.hpp"
#include "OrbitalMechanics.hpp"

using namespace Eigen;
using namespace std;

int main()
{
    // Example input (Keplerian elements)
    double sma = 6852.6;  // Semi-major axis in Km
    double ecc = 0.0013;  // Eccentricity
    double inc = 1.7047;  // Inclination in radians
    double w = 4.7190;    // Argument of Perigee in radians
    double nu = 1.7497;   // True Anomaly in radians
    double raan = 2.7630; // Right Ascension of Ascending Node in radians

    // Get the position and velocity vectors
    auto [r, v] = OrbitalMechanics::keplerian2ijk(sma, ecc, inc, w, nu, raan);

    // Output the results
    std::cout << "Position: [" << r[0] << ", " << r[1] << ", " << r[2] << "] Km" << std::endl;
    std::cout << "Velocity: [" << v[0] << ", " << v[1] << ", " << v[2] << "] Km/s" << std::endl;

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

        cout << "Satellite Number: " << tle.getSatelliteNumber() << endl;
        cout << "Epoch Year: " << tle.getEpochYear() << endl;
        cout << "Mean Motion: " << tle.getMeanMotion() << endl;

        MatrixXd A(2, 3);
        A << 1, 2, 3,
            4, 5, 6;

        MatrixXd B(3, 2);
        B << 7, 8,
            9, 10,
            11, 12;

        MatrixXd result = A * B;
        MatrixXd transposedResult = result.transpose();

        cout << "Result of (A * B)^T:\n"
             << transposedResult << endl;
    }
    catch (const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
