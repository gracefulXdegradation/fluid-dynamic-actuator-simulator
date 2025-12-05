#include <iostream>
#include <cstdlib>
#include <signal.h>
#include "JobProcessor.hpp"
#include "Logger.hpp"

JobProcessor* g_job_processor = nullptr;

void signalHandler(int signal) {
    if (g_job_processor) {
        Logger::info("\nReceived signal " + std::to_string(signal) + ", stopping job processor...");
        g_job_processor->stop();
    }
}

int main(int argc, char* argv[]) {
    // Get database and Redis URLs from environment variables
    const char* db_url = std::getenv("DATABASE_URL");
    const char* redis_url = std::getenv("REDIS_URL");
    
    if (!db_url) {
        Logger::error("Error: DATABASE_URL environment variable not set");
        return 1;
    }
    
    if (!redis_url) {
        Logger::error("Error: REDIS_URL environment variable not set");
        return 1;
    }
    
    // Parse poll interval from command line (optional)
    int poll_interval = 5; // default 5 seconds
    if (argc > 1) {
        try {
            poll_interval = std::stoi(argv[1]);
            if (poll_interval < 1) {
                Logger::error("Warning: Poll interval must be >= 1, using default 5 seconds");
                poll_interval = 5;
            }
        } catch (...) {
            Logger::error("Warning: Invalid poll interval, using default 5 seconds");
        }
    }
    
    try {
        JobProcessor processor(db_url, redis_url);
        g_job_processor = &processor;
        
        // Set up signal handlers for graceful shutdown
        signal(SIGINT, signalHandler);
        signal(SIGTERM, signalHandler);
        
        // Start processing jobs (blocks until stopped)
        processor.start(poll_interval);
        
    } catch (const std::exception& e) {
        Logger::error("Error: " + std::string(e.what()));
        return 1;
    }
    
    return 0;
}

