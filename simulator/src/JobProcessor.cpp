#include "JobProcessor.hpp"
#include "Database.hpp"
#include "RedisClient.hpp"
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
        std::cerr << "JobProcessor is already running" << std::endl;
        return;
    }
    
    running_ = true;
    should_stop_ = false;
    
    std::cout << "JobProcessor started. Polling for jobs every " << poll_interval_seconds << " seconds..." << std::endl;
    std::cout << "Press Ctrl+C to stop." << std::endl;
    
    while (!should_stop_) {
        try {
            // Try to get a job (non-blocking)
            std::string simulation_id = redis_->tryPopJob();
            
            if (!simulation_id.empty()) {
                std::cout << "Found job: " << simulation_id << std::endl;
                processJob(simulation_id);
            } else {
                // No job available, wait before polling again
                std::this_thread::sleep_for(std::chrono::seconds(poll_interval_seconds));
            }
        } catch (const std::exception& e) {
            std::cerr << "Error in job processing loop: " << e.what() << std::endl;
            std::this_thread::sleep_for(std::chrono::seconds(poll_interval_seconds));
        }
    }
    
    running_ = false;
    std::cout << "JobProcessor stopped." << std::endl;
}

void JobProcessor::stop() {
    should_stop_ = true;
    // Wait for the processing loop to finish
    while (running_) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
}

void JobProcessor::processJob(const std::string& simulation_id) {
    std::cout << "Processing simulation: " << simulation_id << std::endl;
    
    // Get the path to the simulator executable
    // Assuming it's in the same directory or we can find it via PATH
    const char* simulator_path = "./build/bin/fds";
    
    // Check if executable exists
    if (access(simulator_path, X_OK) != 0) {
        std::cerr << "Error: Simulator executable not found at " << simulator_path << std::endl;
        db_->updateSimulationStatus(simulation_id, "failed", "Simulator executable not found");
        return;
    }
    
    // Fork and execute the simulator
    pid_t pid = fork();
    
    if (pid < 0) {
        std::cerr << "Error: Failed to fork process" << std::endl;
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
        std::cerr << "Error: Failed to execute simulator" << std::endl;
        exit(1);
    } else {
        // Parent process: wait for child to complete
        int status;
        waitpid(pid, &status, 0);
        
        if (WIFEXITED(status)) {
            int exit_code = WEXITSTATUS(status);
            if (exit_code == 0) {
                std::cout << "Simulation " << simulation_id << " completed successfully" << std::endl;
                // Status should already be updated to "completed" by the simulator
            } else {
                std::cerr << "Simulation " << simulation_id << " failed with exit code " << exit_code << std::endl;
                // Status should already be updated to "failed" by the simulator
            }
        } else if (WIFSIGNALED(status)) {
            int signal = WTERMSIG(status);
            std::cerr << "Simulation " << simulation_id << " was terminated by signal " << signal << std::endl;
            db_->updateSimulationStatus(simulation_id, "failed", "Process terminated by signal");
        }
    }
}

