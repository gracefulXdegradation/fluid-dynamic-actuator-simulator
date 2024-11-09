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

  inline std::chrono::system_clock::time_point parseDateTime(const std::string &datetime_str)
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

  inline std::vector<std::chrono::system_clock::time_point> generateTimePoints(
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

  inline std::string formatTime(std::chrono::system_clock::time_point tp)
  {
    std::time_t time = std::chrono::system_clock::to_time_t(tp);
    std::tm *tm = std::localtime(&time);

    std::ostringstream oss;
    oss << std::put_time(tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
  }

} // namespace DateTime

#endif // DATETIME_HPP
