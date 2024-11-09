#include <iostream>
#include <string>
#include <chrono>
#include <ctime>
#include <vector>
#include <iomanip>
#include <Eigen/Dense>
#include "TLE.h"
#include "config.hpp"

using namespace Eigen;
using namespace std;

// Helper function to format the time_point as a string for printing
std::string format_time(const std::chrono::system_clock::time_point &time_point)
{
    std::time_t time = std::chrono::system_clock::to_time_t(time_point);
    std::tm tm = *std::localtime(&time);

    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", &tm);
    return std::string(buffer);
}

// Function to generate discrete time points
std::vector<std::chrono::system_clock::time_point> generateTimePoints(
    const std::chrono::system_clock::time_point &start_date_time,
    const std::chrono::system_clock::time_point &end_date_time,
    const std::chrono::milliseconds &control_time_step)
{
    std::vector<std::chrono::system_clock::time_point> date_times;
    for (auto current = start_date_time; current <= end_date_time; current += control_time_step)
    {
        date_times.push_back(current);
    }
    return date_times;
}

int main()
{
    // Create an instance of ConfigParser with the path to your config.json
    Config config(std::string(BUILD_OUTPUT_PATH) + "/config.json");

    // Generate discrete time points using the function
    auto date_times = generateTimePoints(config.getStartDateTime(), config.getEndDateTime(), config.getControlTimeStep());

    // Output the time points
    std::cout << "Executing simulation" << std::endl;
    std::cout << "====================" << std::endl;
    std::cout << "Start date: " << format_time(date_times.at(0)) << std::endl;

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
