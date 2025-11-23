#ifndef REDIS_CLIENT_HPP
#define REDIS_CLIENT_HPP

#include <string>
#include <memory>

// Forward declaration
struct redisContext;

class RedisClient {
public:
    RedisClient(const std::string& redis_url);
    ~RedisClient();

    // Test connection
    bool testConnection();

    // Push a simulation job to the queue
    bool pushJob(const std::string& simulation_id);

    // Blocking pop from the queue (waits for a job)
    // Returns empty string if connection lost or error
    std::string popJob(int timeout_seconds = 0);

    // Non-blocking pop from the queue
    // Returns empty string if no job available
    std::string tryPopJob();

    // Get queue length
    int getQueueLength();

private:
    std::unique_ptr<redisContext> context_;
    std::string queue_name_;
    std::string redis_url_;

    // Parse Redis URL (redis://host:port or redis://host:port/db)
    void parseRedisUrl(const std::string& url, std::string& host, int& port, int& db);
};

#endif // REDIS_CLIENT_HPP

