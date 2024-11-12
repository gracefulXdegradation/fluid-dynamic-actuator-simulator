#include <chrono>
#include <fstream>
#include <Eigen/Dense>
#include "json.hpp"
#include "dateTime.hpp"

using json = nlohmann::json;
using namespace std::chrono;

class Config
{
public:
  Config(const std::string &config_file)
  {
    std::ifstream file(config_file);
    if (!file.is_open())
    {
      throw std::runtime_error("Failed to open the config file");
    }

    file >> config;

    start_date_time = DateTime::parseDateTime(config["start_date_time"]);
    end_date_time = DateTime::parseDateTime(config["end_date_time"]);
    control_time_step = milliseconds(config["control_time_step"]);

    if (config.contains("ground_station_lla") && config["ground_station_lla"].is_object())
    {
      double lat = config["ground_station_lla"]["lat"];
      double lon = config["ground_station_lla"]["lon"];
      double alt = config["ground_station_lla"]["alt"];
      gs_r.resize(3);
      gs_r << lat, lon, alt;
    }
    else
    {
      throw std::runtime_error("Field 'ground_station_lla' not found or is not an object.");
    }
  }

  system_clock::time_point getStartDateTime() const
  {
    return start_date_time;
  }

  system_clock::time_point getEndDateTime() const
  {
    return end_date_time;
  }

  milliseconds getControlTimeStep() const
  {
    return control_time_step;
  }

  const Eigen::VectorXd &getGroundStationPosition() const
  {
    return gs_r;
  }

private:
  json config;
  system_clock::time_point start_date_time;
  system_clock::time_point end_date_time;
  milliseconds control_time_step;
  Eigen::VectorXd gs_r;
};