#include <iostream>
#include <fstream>
#include <string>

#include "TLE.h"

using namespace std;

int main()
{
    try
    {
        TLE tle = TLE::fromFile(std::string(BUILD_OUTPUT_PATH) + "/tle.txt");

        std::cout << "Satellite Number: " << tle.getSatelliteNumber() << std::endl;
        std::cout << "Epoch Year: " << tle.getEpochYear() << std::endl;
        std::cout << "Mean Motion: " << tle.getMeanMotion() << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}
