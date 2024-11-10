#ifndef DATETIME_HPP
#define DATETIME_HPP

#include <chrono>
#include <string>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <stdexcept>
#include <vector>

namespace DateTime
{
  using namespace std::chrono;

  constexpr std::chrono::seconds SECONDS_IN_A_DAY(86400);

  inline system_clock::time_point parseDateTime(const std::string &datetime_str)
  {
    std::tm tm = {};
    std::istringstream ss(datetime_str);
    ss >> std::get_time(&tm, "%Y-%m-%d %H:%M:%S");

    if (ss.fail())
    {
      throw std::runtime_error("Failed to parse date-time string");
    }

    return system_clock::from_time_t(std::mktime(&tm));
  }

  inline std::vector<system_clock::time_point> generateTimePoints(
      const system_clock::time_point &start_date_time,
      const system_clock::time_point &end_date_time,
      const milliseconds &control_time_step)
  {
    std::vector<system_clock::time_point> date_times;
    for (auto current = start_date_time; current <= end_date_time; current += control_time_step)
    {
      date_times.push_back(current);
    }
    return date_times;
  }

  inline std::string formatTime(system_clock::time_point tp)
  {
    std::time_t time = system_clock::to_time_t(tp);
    std::tm *tm = std::localtime(&time);

    std::ostringstream oss;
    oss << std::put_time(tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
  }

  inline std::string getCurrentTimestamp()
  {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y%m%d_%H%M%S");
    return ss.str();
  }

}

#endif
