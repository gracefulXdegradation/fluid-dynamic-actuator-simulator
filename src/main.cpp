#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <ctime>
#include <vector>
#include <iomanip>
#include <Eigen/Dense>
#include "TLE.h"
#include "json.hpp"

using json = nlohmann::json;
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

// Function to parse date-time string to time_point
std::chrono::system_clock::time_point parseDateTime(const std::string &datetime_str)
{
    std::tm tm = {};
    std::istringstream ss(datetime_str);
    ss >> std::get_time(&tm, "%Y-%m-%d %H:%M:%S");
    if (ss.fail())
    {
        throw std::runtime_error("Failed to parse date-time string");
    }
    return std::chrono::system_clock::from_time_t(std::mktime(&tm));
}

int main()
{
    using namespace std::chrono;

    // Read JSON configuration file
    std::ifstream file(std::string(BUILD_OUTPUT_PATH) + "/config.json");
    json config;
    file >> config;

    // Parse start and end date-times from JSON
    auto start_date_time = parseDateTime(config["start_date_time"].get<std::string>());
    auto end_date_time = parseDateTime(config["end_date_time"].get<std::string>());

    // Define control time step of 0.5 seconds (500 milliseconds)
    auto control_time_step = milliseconds(500);

    // Generate discrete time points using the function
    auto date_times = generateTimePoints(start_date_time, end_date_time, control_time_step);

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
