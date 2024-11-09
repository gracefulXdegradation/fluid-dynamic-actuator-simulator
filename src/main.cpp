#include <iostream>
#include <string>
#include <chrono>
#include <ctime>
#include <vector>
#include <iomanip>
#include <Eigen/Dense>
#include "TLE.h"
#include "config.hpp"
#include "dateTime.hpp"

using namespace Eigen;
using namespace std;

int main()
{
    // Create an instance of ConfigParser with the path to your config.json
    Config config(std::string(BUILD_OUTPUT_PATH) + "/config.json");

    // Generate discrete time points using the function
    auto date_times = DateTime::generateTimePoints(config.getStartDateTime(), config.getEndDateTime(), config.getControlTimeStep());

    // Output the time points
    std::cout << "Executing simulation" << std::endl;
    std::cout << "====================" << std::endl;
    std::cout << "Start date: " << DateTime::formatTime(date_times.at(0)) << std::endl;

    try
    {
        TLE tle = TLE::fromFile(std::string(BUILD_OUTPUT_PATH) + "/tle.txt");

        std::cout << "Satellite Number: " << tle.getSatelliteNumber() << std::endl;
        std::cout << "Epoch Year: " << tle.getEpochYear() << std::endl;
        std::cout << "Mean Motion: " << tle.getMeanMotion() << std::endl;

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
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
