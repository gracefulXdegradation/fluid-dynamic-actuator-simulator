#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <chrono>

class Logger {
public:
    // Get current timestamp as string in format: YYYY-MM-DD HH:MM:SS.mmm
    static std::string getCurrentTimestamp() {
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()) % 1000;
        
        std::stringstream ss;
        ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
        ss << '.' << std::setfill('0') << std::setw(3) << ms.count();
        return ss.str();
    }

    // Log info message with timestamp to stdout
    static void info(const std::string& message) {
        std::cout << "[" << getCurrentTimestamp() << "] " << message << std::endl;
    }

    // Log error message with timestamp to stderr
    static void error(const std::string& message) {
        std::cerr << "[" << getCurrentTimestamp() << "] " << message << std::endl;
    }
};

#endif // LOGGER_HPP

