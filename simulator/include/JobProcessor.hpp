#ifndef JOB_PROCESSOR_HPP
#define JOB_PROCESSOR_HPP

#include <string>
#include <atomic>
#include <memory>

class Database;
class RedisClient;

class JobProcessor {
public:
    JobProcessor(const std::string& database_url, const std::string& redis_url);
    ~JobProcessor();

    // Start processing jobs (blocks until stopped)
    void start(int poll_interval_seconds = 5);

    // Stop processing jobs
    void stop();

    // Check if processor is running
    bool isRunning() const { return running_; }

private:
    std::unique_ptr<Database> db_;
    std::unique_ptr<RedisClient> redis_;
    std::atomic<bool> running_;
    std::atomic<bool> should_stop_;

    // Process a single simulation job
    void processJob(const std::string& simulation_id);
};

#endif // JOB_PROCESSOR_HPP

