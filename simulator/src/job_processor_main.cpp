#include <iostream>
#include <cstdlib>
#include <signal.h>
#include "JobProcessor.hpp"

JobProcessor* g_job_processor = nullptr;

void signalHandler(int signal) {
    if (g_job_processor) {
        std::cout << "\nReceived signal " << signal << ", stopping job processor..." << std::endl;
        g_job_processor->stop();
    }
}

int main(int argc, char* argv[]) {
    // Get database and Redis URLs from environment variables
    const char* db_url = std::getenv("DATABASE_URL");
    const char* redis_url = std::getenv("REDIS_URL");
    
    if (!db_url) {
        std::cerr << "Error: DATABASE_URL environment variable not set" << std::endl;
        return 1;
    }
    
    if (!redis_url) {
        std::cerr << "Error: REDIS_URL environment variable not set" << std::endl;
        return 1;
    }
    
    // Parse poll interval from command line (optional)
    int poll_interval = 5; // default 5 seconds
    if (argc > 1) {
        try {
            poll_interval = std::stoi(argv[1]);
            if (poll_interval < 1) {
                std::cerr << "Warning: Poll interval must be >= 1, using default 5 seconds" << std::endl;
                poll_interval = 5;
            }
        } catch (...) {
            std::cerr << "Warning: Invalid poll interval, using default 5 seconds" << std::endl;
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
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

