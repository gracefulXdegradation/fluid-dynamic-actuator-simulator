#include "JobProcessor.hpp"
#include "Database.hpp"
#include "RedisClient.hpp"
#include "Logger.hpp"
#include <iostream>
#include <thread>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <sys/wait.h>

JobProcessor::JobProcessor(const std::string& database_url, const std::string& redis_url)
    : running_(false), should_stop_(false) {
    try {
        db_ = std::make_unique<Database>(database_url);
        if (!db_->testConnection()) {
            throw std::runtime_error("Failed to connect to database");
        }
        
        redis_ = std::make_unique<RedisClient>(redis_url);
        if (!redis_->testConnection()) {
            throw std::runtime_error("Failed to connect to Redis");
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("JobProcessor initialization failed: " + std::string(e.what()));
    }
}

JobProcessor::~JobProcessor() {
    stop();
}

void JobProcessor::start(int poll_interval_seconds) {
    if (running_) {
        Logger::error("JobProcessor is already running");
        return;
    }
    
    running_ = true;
    should_stop_ = false;
    
    Logger::info("JobProcessor started. Polling for jobs every " + std::to_string(poll_interval_seconds) + " seconds...");
    Logger::info("Press Ctrl+C to stop.");
    
    while (!should_stop_) {
        try {
            // Try to get a job (non-blocking)
            std::string simulation_id = redis_->tryPopJob();
            
            if (!simulation_id.empty()) {
                Logger::info("\nFound job: " + simulation_id);
                processJob(simulation_id);
            } else {
                // No job available, wait before polling again
                std::this_thread::sleep_for(std::chrono::seconds(poll_interval_seconds));
            }
        } catch (const std::exception& e) {
            Logger::error("Error in job processing loop: " + std::string(e.what()));
            std::this_thread::sleep_for(std::chrono::seconds(poll_interval_seconds));
        }
    }
    
    running_ = false;
    Logger::info("JobProcessor stopped.");
}

void JobProcessor::stop() {
    should_stop_ = true;
    // Wait for the processing loop to finish
    while (running_) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
}

void JobProcessor::processJob(const std::string& simulation_id) {
    Logger::info("Processing simulation: " + simulation_id);
    
    // Get the path to the simulator executable
    // Assuming it's in the same directory or we can find it via PATH
    const char* simulator_path = "./build/bin/fds";
    
    // Check if executable exists
    if (access(simulator_path, X_OK) != 0) {
        Logger::error("Error: Simulator executable not found at " + std::string(simulator_path));
        db_->updateSimulationStatus(simulation_id, "failed", "Simulator executable not found");
        return;
    }
    
    // Fork and execute the simulator
    pid_t pid = fork();
    
    if (pid < 0) {
        Logger::error("Error: Failed to fork process");
        db_->updateSimulationStatus(simulation_id, "failed", "Failed to fork process");
        return;
    }
    
    if (pid == 0) {
        // Child process: execute the simulator
        char* args[] = {
            const_cast<char*>(simulator_path),
            const_cast<char*>(simulation_id.c_str()),
            nullptr
        };
        
        execvp(simulator_path, args);
        
        // If execvp returns, it failed
        Logger::error("Error: Failed to execute simulator");
        exit(1);
    } else {
        // Parent process: wait for child to complete
        int status;
        waitpid(pid, &status, 0);
        
        if (WIFEXITED(status)) {
            int exit_code = WEXITSTATUS(status);
            if (exit_code == 0) {
                Logger::info("Simulation " + simulation_id + " completed successfully");
                // Status should already be updated to "completed" by the simulator
            } else {
                Logger::error("Simulation " + simulation_id + " failed with exit code " + std::to_string(exit_code));
                // Status should already be updated to "failed" by the simulator
            }
        } else if (WIFSIGNALED(status)) {
            int signal = WTERMSIG(status);
            Logger::error("Simulation " + simulation_id + " was terminated by signal " + std::to_string(signal));
            db_->updateSimulationStatus(simulation_id, "failed", "Process terminated by signal");
        }
    }
}

