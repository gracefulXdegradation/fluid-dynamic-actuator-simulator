#include "TLE.h"
#include <fstream>
#include <stdexcept>

TLE::TLE(const std::string &line1, const std::string &line2)
{
  parseTLE(line1, line2);
}

void TLE::parseTLE(const std::string &line1, const std::string &line2)
{
  // Parse Line 1
  lineNumber1 = std::stoi(line1.substr(0, 1));
  satelliteNumber1 = std::stoi(line1.substr(2, 5));
  classification = line1[7];
  launchYear = std::stoi(line1.substr(9, 2));
  launchNumber = std::stoi(line1.substr(11, 3));
  launchPiece = line1.substr(14, 3);
  epochYear = std::stoi(line1.substr(18, 2));
  epochDay = std::stod(line1.substr(20, 12));
  firstDerivMeanMotion = std::stod(line1.substr(33, 10));
  secondDerivMeanMotion = std::stod(line1.substr(44, 8)) * 1e-5;
  bstarDrag = std::stod(line1.substr(53, 8)) * 1e-5;
  ephemerisType = std::stoi(line1.substr(62, 1));
  elementSetNumber = std::stoi(line1.substr(64, 4));
  checksum1 = std::stoi(line1.substr(68, 1));

  // Parse Line 2
  lineNumber2 = std::stoi(line2.substr(0, 1));
  satelliteNumber2 = std::stoi(line2.substr(2, 5));
  inclination = std::stod(line2.substr(8, 8));
  rightAscension = std::stod(line2.substr(17, 8));
  eccentricity = std::stod("0." + line2.substr(26, 7));
  argumentOfPerigee = std::stod(line2.substr(34, 8));
  meanAnomaly = std::stod(line2.substr(43, 8));
  meanMotion = std::stod(line2.substr(52, 11));
  revolutionNumber = std::stoi(line2.substr(63, 5));
  checksum2 = std::stoi(line2.substr(68, 1));
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
