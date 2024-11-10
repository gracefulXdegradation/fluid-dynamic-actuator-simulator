#ifndef TLE_H
#define TLE_H

#include <string>
#include <chrono>

class TLE
{
public:
  // Constructor that accepts two lines of TLE data
  TLE(const std::string &line1, const std::string &line2);

  // Static function to create a TLE object from a file
  static TLE fromFile(const std::string &filePath);

  std::chrono::system_clock::time_point getEpoch() const { return epoch; }
  double getMeanAnomaly() const { return meanAnomaly; }
  double getMeanMotion() const { return meanMotion; }
  double getRightAscension() const { return rightAscension; }
  double getEccentricity() const { return eccentricity; }
  double getInclination() const { return inclination; }
  double getArgumentOfPerigee() const { return argumentOfPerigee; }
  double getSemiMajorAxis() const { return semiMajorAxis; }

  static const double GM; // gravitational parameter

private:
  // Line 1 fields
  int noradId;
  char classification;
  int launchYear;
  int launchNumber;
  std::string launchPiece;
  std::chrono::system_clock::time_point epoch;
  double firstDerivMeanMotion;
  double secondDerivMeanMotion;
  double bstarDrag;
  int ephemerisType;
  int elementSetNumber;
  int checksum1;

  // Line 2 fields
  int satelliteNumber2;
  double inclination;
  double rightAscension;
  double eccentricity;
  double argumentOfPerigee;
  double meanAnomaly;
  double meanMotion;
  int revolutionNumber;
  int checksum2;

  double semiMajorAxis; // km

  void parseTLE(const std::string &line1, const std::string &line2);
};

#endif
