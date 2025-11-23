#ifndef DATABASE_HPP
#define DATABASE_HPP

#include <string>
#include <chrono>
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "json.hpp"
#include <pqxx/pqxx>

using json = nlohmann::json;
using namespace std::chrono;

struct SimulationParams {
    std::string tle_line1;
    std::string tle_line2;
    system_clock::time_point start_date_time;
    system_clock::time_point end_date_time;
    milliseconds control_time_step;
    Eigen::Vector3d ground_station_lla;
    double ground_station_elevation;
};

struct SimulationMetric {
    int step_index;
    int64_t timestamp; // milliseconds since epoch
    Eigen::Vector3d euler_angles;
    Eigen::Vector3d ang_mom_body_frame;
    Eigen::Vector4d a_control_torque;
    Eigen::Vector4d a_command;
    Eigen::Matrix<double, 15, 1> state;
    double distance;
};

class Database {
public:
    Database(const std::string& connection_string);
    ~Database();

    // Read simulation parameters from database
    SimulationParams getSimulationParams(const std::string& simulation_id);

    // Update simulation status
    void updateSimulationStatus(const std::string& simulation_id, 
                                const std::string& status,
                                const std::string& error_message = "");

    // Write metrics to database (batch insert for efficiency)
    void writeMetrics(const std::string& simulation_id, 
                      const std::vector<SimulationMetric>& metrics);

    // Create or update simulation parameters in database
    void createOrUpdateSimulationParams(const std::string& simulation_id, 
                                        const json& params_json);

    // Test database connection
    bool testConnection();

private:
    std::unique_ptr<pqxx::connection> conn_;
    std::string connection_string_;

    // Helper to parse JSON input_parameters
    SimulationParams parseInputParameters(const json& params_json);
};

#endif // DATABASE_HPP

