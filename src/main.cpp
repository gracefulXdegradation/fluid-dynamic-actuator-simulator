#include <iostream>
#include <fstream>
#include <string>

#include "Hello.h"

using namespace std;

struct TLE
{
    // Line 1 fields
    int lineNumber1;
    int satelliteNumber1;
    char classification;
    int launchYear;
    int launchNumber;
    std::string launchPiece;
    int epochYear;
    double epochDay;
    double firstDerivMeanMotion;
    double secondDerivMeanMotion;
    double bstarDrag;
    int ephemerisType;
    int elementSetNumber;
    int checksum1;

    // Line 2 fields
    int lineNumber2;
    int satelliteNumber2;
    double inclination;
    double rightAscension;
    double eccentricity;
    double argumentOfPerigee;
    double meanAnomaly;
    double meanMotion;
    int revolutionNumber;
    int checksum2;
};

TLE parseTLE(const std::string &line1, const std::string &line2)
{
    TLE tle;

    // Parse Line 1
    tle.lineNumber1 = std::stoi(line1.substr(0, 1));
    tle.satelliteNumber1 = std::stoi(line1.substr(2, 5));
    tle.classification = line1[7];
    tle.launchYear = std::stoi(line1.substr(9, 2));
    tle.launchNumber = std::stoi(line1.substr(11, 3));
    tle.launchPiece = line1.substr(14, 3);
    tle.epochYear = std::stoi(line1.substr(18, 2));
    tle.epochDay = std::stod(line1.substr(20, 12));
    tle.firstDerivMeanMotion = std::stod(line1.substr(33, 10));
    tle.secondDerivMeanMotion = std::stod(line1.substr(44, 8)) * 1e-5;
    tle.bstarDrag = std::stod(line1.substr(53, 8)) * 1e-5;
    tle.ephemerisType = std::stoi(line1.substr(62, 1));
    tle.elementSetNumber = std::stoi(line1.substr(64, 4));
    tle.checksum1 = std::stoi(line1.substr(68, 1));

    // Parse Line 2
    tle.lineNumber2 = std::stoi(line2.substr(0, 1));
    tle.satelliteNumber2 = std::stoi(line2.substr(2, 5));
    tle.inclination = std::stod(line2.substr(8, 8));
    tle.rightAscension = std::stod(line2.substr(17, 8));
    tle.eccentricity = std::stod("0." + line2.substr(26, 7)); // Assumes decimal point before
    tle.argumentOfPerigee = std::stod(line2.substr(34, 8));
    tle.meanAnomaly = std::stod(line2.substr(43, 8));
    tle.meanMotion = std::stod(line2.substr(52, 11));
    tle.revolutionNumber = std::stoi(line2.substr(63, 5));
    tle.checksum2 = std::stoi(line2.substr(68, 1));

    return tle;
}

int main()
{
    std::ifstream file(std::string(BUILD_OUTPUT_PATH) + "/tle.txt");
    if (!file.is_open())
    {
        std::cerr << "Failed to open file." << std::endl;
        return 1;
    }

    std::string line1, line2;
    if (std::getline(file, line1) && std::getline(file, line2))
    {
        TLE tle = parseTLE(line1, line2);

        // Output parsed data for verification
        std::cout << "Parsed TLE Data:\n";
        std::cout << "Satellite Number: " << tle.satelliteNumber1 << "\n";
        std::cout << "Epoch Year: " << tle.epochYear << "\n";
        std::cout << "Epoch Day: " << tle.epochDay << "\n";
        std::cout << "Inclination: " << tle.inclination << "\n";
        std::cout << "Right Ascension: " << tle.rightAscension << "\n";
        std::cout << "Eccentricity: " << tle.eccentricity << "\n";
        std::cout << "Mean Motion: " << tle.meanMotion << "\n";
        std::cout << "Revolution Number: " << tle.revolutionNumber << "\n";
    }
    else
    {
        std::cerr << "File does not contain enough lines for a TLE." << std::endl;
    }

    file.close();

    Hello hi;
    hi.print();
    return 0;
}
