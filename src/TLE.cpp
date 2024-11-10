#include "TLE.h"
#include "MathHelpers.hpp"
#include "dateTime.hpp"
#include <fstream>
#include <stdexcept>
#include <chrono>
#include <ctime>

using namespace MathHelpers;
using namespace std::chrono;

TLE::TLE(const std::string &line1, const std::string &line2)
{
  parseTLE(line1, line2);
}

int getFullYear(int yearTwoDigits, int cutoff = 57)
{
  if (yearTwoDigits >= 0 && yearTwoDigits <= 99)
  {
    return (yearTwoDigits >= cutoff) ? (1900 + yearTwoDigits) : (2000 + yearTwoDigits);
  }
  else
  {
    throw std::invalid_argument("Year must be a two-digit number between 0 and 99.");
  }
}

system_clock::time_point dayOfYearToDateTime(double dayOfYear, int year)
{
  std::tm startOfTheYear = {0, 0, 0, 1, 0, year - 1900};
  system_clock::time_point yearStart = system_clock::from_time_t(std::mktime(&startOfTheYear));
  duration<int> secondsFromYearStartTilTheMoment(static_cast<int>(dayOfYear * DateTime::SECONDS_IN_A_DAY.count()));
  return yearStart + secondsFromYearStartTilTheMoment;
}

void TLE::parseTLE(const std::string &line1, const std::string &line2)
{
  // Parse Line 1
  noradId = stoi(line1.substr(2, 5));
  classification = line1[7];
  launchYear = stoi(line1.substr(9, 2));
  launchNumber = stoi(line1.substr(11, 3));
  launchPiece = line1.substr(14, 3);
  epoch = dayOfYearToDateTime(stod(line1.substr(20, 12)), getFullYear(stoi(line1.substr(18, 2))));
  firstDerivMeanMotion = stod(line1.substr(33, 10));
  secondDerivMeanMotion = stod(line1.substr(44, 8)) * 1e-5;
  bstarDrag = stod(line1.substr(53, 8)) * 1e-5;
  ephemerisType = stoi(line1.substr(62, 1));
  elementSetNumber = stoi(line1.substr(64, 4));
  checksum1 = stoi(line1.substr(68, 1));

  // Parse Line 2
  satelliteNumber2 = stoi(line2.substr(2, 5));
  inclination = deg2rad(stod(line2.substr(8, 8)));
  rightAscension = deg2rad(stod(line2.substr(17, 8)));
  eccentricity = stod("0." + line2.substr(26, 7));
  argumentOfPerigee = deg2rad(stod(line2.substr(34, 8)));
  meanAnomaly = deg2rad(stod(line2.substr(43, 8)));
  meanMotion = 2 * M_PI * stod(line2.substr(52, 11)) / DateTime::SECONDS_IN_A_DAY.count();
  revolutionNumber = stoi(line2.substr(63, 5));
  checksum2 = stoi(line2.substr(68, 1));
}

TLE TLE::fromFile(const std::string &filePath)
{
  std::ifstream file(filePath);
  if (!file)
  {
    throw std::runtime_error("Unable to open TLE file: " + filePath);
  }

  std::string line1, line2;

  // Read two lines (one TLE set)
  if (!std::getline(file, line1) || !std::getline(file, line2))
  {
    throw std::runtime_error("File does not contain enough lines for a TLE set");
  }

  return TLE(line1, line2);
}
