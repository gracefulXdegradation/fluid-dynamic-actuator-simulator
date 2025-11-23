#include "Database.hpp"
#include "dateTime.hpp"
#include <pqxx/pqxx>
#include <stdexcept>
#include <sstream>
#include <iostream>

Database::Database(const std::string& connection_string) 
    : connection_string_(connection_string) {
    try {
        conn_ = std::make_unique<pqxx::connection>(connection_string);
        if (!conn_->is_open()) {
            throw std::runtime_error("Failed to open database connection");
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("Database connection error: " + std::string(e.what()));
    }
}

Database::~Database() {
    // Connection will be automatically closed when unique_ptr is destroyed
}

bool Database::testConnection() {
    try {
        if (!conn_ || !conn_->is_open()) {
            return false;
        }
        pqxx::work txn(*conn_);
        txn.exec("SELECT 1");
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Database connection test failed: " << e.what() << std::endl;
        return false;
    }
}

SimulationParams Database::getSimulationParams(const std::string& simulation_id) {
    try {
        pqxx::work txn(*conn_);
        
        // Query simulation parameters
        std::string query = "SELECT input_parameters FROM simulations WHERE id = $1";
        pqxx::result result = txn.exec_params(query, simulation_id);
        
        if (result.empty()) {
            throw std::runtime_error("Simulation with id " + simulation_id + " not found");
        }
        
        // Get the JSONB field
        std::string params_json_str = result[0][0].as<std::string>();
        json params_json = json::parse(params_json_str);
        
        return parseInputParameters(params_json);
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to get simulation parameters: " + std::string(e.what()));
    }
}

SimulationParams Database::parseInputParameters(const json& params_json) {
    SimulationParams params;
    
    // Parse TLE lines
    if (!params_json.contains("tle_line1") || !params_json.contains("tle_line2")) {
        throw std::runtime_error("TLE lines not found in input_parameters");
    }
    params.tle_line1 = params_json["tle_line1"].get<std::string>();
    params.tle_line2 = params_json["tle_line2"].get<std::string>();
    
    // Parse date times
    if (!params_json.contains("start_date_time") || !params_json.contains("end_date_time")) {
        throw std::runtime_error("Date times not found in input_parameters");
    }
    params.start_date_time = DateTime::parseDateTime(params_json["start_date_time"].get<std::string>());
    params.end_date_time = DateTime::parseDateTime(params_json["end_date_time"].get<std::string>());
    
    // Parse control time step
    if (!params_json.contains("control_time_step")) {
        throw std::runtime_error("control_time_step not found in input_parameters");
    }
    params.control_time_step = milliseconds(params_json["control_time_step"].get<int>());
    
    // Parse ground station parameters
    if (!params_json.contains("ground_station_lla") || !params_json["ground_station_lla"].is_object()) {
        throw std::runtime_error("ground_station_lla not found or is not an object");
    }
    params.ground_station_lla << 
        params_json["ground_station_lla"]["lat"].get<double>(),
        params_json["ground_station_lla"]["lon"].get<double>(),
        params_json["ground_station_lla"]["alt"].get<double>();
    
    if (!params_json.contains("ground_station_elevation")) {
        throw std::runtime_error("ground_station_elevation not found in input_parameters");
    }
    params.ground_station_elevation = params_json["ground_station_elevation"].get<double>();
    
    return params;
}

void Database::updateSimulationStatus(const std::string& simulation_id, 
                                      const std::string& status,
                                      const std::string& error_message) {
    try {
        pqxx::work txn(*conn_);
        
        std::string query;
        if (status == "running") {
            query = "UPDATE simulations SET status = $1, started_at = CURRENT_TIMESTAMP WHERE id = $2";
            txn.exec_params(query, status, simulation_id);
        } else if (status == "completed") {
            query = "UPDATE simulations SET status = $1, completed_at = CURRENT_TIMESTAMP WHERE id = $2";
            txn.exec_params(query, status, simulation_id);
        } else if (status == "failed") {
            query = "UPDATE simulations SET status = $1, completed_at = CURRENT_TIMESTAMP, error_message = $3 WHERE id = $2";
            txn.exec_params(query, status, simulation_id, error_message);
        } else {
            query = "UPDATE simulations SET status = $1 WHERE id = $2";
            txn.exec_params(query, status, simulation_id);
        }
        
        txn.commit();
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to update simulation status: " + std::string(e.what()));
    }
}

void Database::writeMetrics(const std::string& simulation_id, 
                           const std::vector<SimulationMetric>& metrics) {
    if (metrics.empty()) {
        return;
    }
    
    try {
        pqxx::work txn(*conn_);
        
        // Prepare batch insert using prepared statement approach
        // We'll use a single transaction with multiple inserts for efficiency
        std::string insert_query = 
            "INSERT INTO simulation_metrics (simulation_id, step_index, timestamp, "
            "euler_angles, ang_mom_body_frame, a_control_torque, a_command, state, distance) "
            "VALUES ($1, $2, $3, $4, $5, $6, $7, $8, $9) "
            "ON CONFLICT (simulation_id, step_index) DO NOTHING";
        
        // Use pipeline for batch inserts (more efficient)
        for (const auto& metric : metrics) {
            // Convert Eigen vectors to PostgreSQL array format
            std::string euler_angles_str = "{" + 
                std::to_string(metric.euler_angles(0)) + "," +
                std::to_string(metric.euler_angles(1)) + "," +
                std::to_string(metric.euler_angles(2)) + "}";
            
            std::string ang_mom_str = "{" + 
                std::to_string(metric.ang_mom_body_frame(0)) + "," +
                std::to_string(metric.ang_mom_body_frame(1)) + "," +
                std::to_string(metric.ang_mom_body_frame(2)) + "}";
            
            std::string torque_str = "{" + 
                std::to_string(metric.a_control_torque(0)) + "," +
                std::to_string(metric.a_control_torque(1)) + "," +
                std::to_string(metric.a_control_torque(2)) + "," +
                std::to_string(metric.a_control_torque(3)) + "}";
            
            std::string command_str = "{" + 
                std::to_string(metric.a_command(0)) + "," +
                std::to_string(metric.a_command(1)) + "," +
                std::to_string(metric.a_command(2)) + "," +
                std::to_string(metric.a_command(3)) + "}";
            
            std::string state_str = "{";
            for (int i = 0; i < 15; i++) {
                state_str += std::to_string(metric.state(i));
                if (i < 14) state_str += ",";
            }
            state_str += "}";
            
            txn.exec_params(insert_query,
                simulation_id,
                metric.step_index,
                metric.timestamp,
                euler_angles_str,
                ang_mom_str,
                torque_str,
                command_str,
                state_str,
                metric.distance
            );
        }
        
        txn.commit();
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to write metrics: " + std::string(e.what()));
    }
}

void Database::createOrUpdateSimulationParams(const std::string& simulation_id, 
                                              const json& params_json) {
    try {
        pqxx::work txn(*conn_);
        
        // Convert JSON to string
        std::string params_json_str = params_json.dump();
        
        // Use INSERT ... ON CONFLICT to create or update
        // Only update input_parameters, leave status unchanged
        std::string query = 
            "INSERT INTO simulations (id, input_parameters) "
            "VALUES ($1, $2::jsonb) "
            "ON CONFLICT (id) DO UPDATE SET input_parameters = $2::jsonb";
        
        txn.exec_params(query, simulation_id, params_json_str);
        txn.commit();
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to create/update simulation parameters: " + std::string(e.what()));
    }
}

