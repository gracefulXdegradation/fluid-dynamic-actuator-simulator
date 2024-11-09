#include <chrono>
#include <fstream>
#include "json.hpp"
#include "dateTime.hpp"

using json = nlohmann::json;

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
    control_time_step = std::chrono::milliseconds(config["control_time_step"]);
  }

  std::chrono::system_clock::time_point getStartDateTime() const
  {
    return start_date_time;
  }

  std::chrono::system_clock::time_point getEndDateTime() const
  {
    return end_date_time;
  }

  std::chrono::milliseconds getControlTimeStep() const
  {
    return control_time_step;
  }

private:
  json config;
  std::chrono::system_clock::time_point start_date_time;
  std::chrono::system_clock::time_point end_date_time;
  std::chrono::milliseconds control_time_step;
};