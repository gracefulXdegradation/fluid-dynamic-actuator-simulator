#ifndef TLE_H
#define TLE_H

#include <string>

class TLE
{
public:
  // Constructor that accepts two lines of TLE data
  TLE(const std::string &line1, const std::string &line2);

  // Static function to create a TLE object from a file
  static TLE fromFile(const std::string &filePath);

  // Getters for accessing fields (optional based on requirements)
  int getSatelliteNumber() const { return satelliteNumber1; }
  double getMeanMotion() const { return meanMotion; }
  int getEpochYear() const { return epochYear; }
  // Add other getters as needed

private:
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

  // Private method to parse the two lines of TLE data
  void parseTLE(const std::string &line1, const std::string &line2);
};

#endif // TLE_H
