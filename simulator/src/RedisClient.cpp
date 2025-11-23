#include "RedisClient.hpp"
#include <hiredis/hiredis.h>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <cstring>

RedisClient::RedisClient(const std::string& redis_url) 
    : redis_url_(redis_url), queue_name_("simulation_jobs") {
    
    std::string host;
    int port;
    int db;
    parseRedisUrl(redis_url, host, port, db);
    
    // Connect to Redis
    struct timeval timeout = { 1, 500000 }; // 1.5 seconds
    redisContext* ctx = redisConnectWithTimeout(host.c_str(), port, timeout);
    
    if (ctx == nullptr || ctx->err) {
        if (ctx) {
            throw std::runtime_error("Redis connection error: " + std::string(ctx->errstr));
        } else {
            throw std::runtime_error("Redis connection error: Can't allocate redis context");
        }
    }
    
    context_ = std::unique_ptr<redisContext>(ctx);
    
    // Select database if specified
    if (db > 0) {
        redisReply* reply = (redisReply*)redisCommand(context_.get(), "SELECT %d", db);
        if (reply == nullptr) {
            throw std::runtime_error("Redis SELECT command failed");
        }
        freeReplyObject(reply);
    }
}

RedisClient::~RedisClient() {
    if (context_) {
        redisFree(context_.get());
    }
}

void RedisClient::parseRedisUrl(const std::string& url, std::string& host, int& port, int& db) {
    // Parse redis://host:port or redis://host:port/db
    // Default: localhost:6379, db 0
    
    host = "localhost";
    port = 6379;
    db = 0;
    
    if (url.empty() || url.find("redis://") != 0) {
        return; // Use defaults
    }
    
    std::string remaining = url.substr(8); // Skip "redis://"
    
    // Find host:port
    size_t colon_pos = remaining.find(':');
    size_t slash_pos = remaining.find('/');
    
    if (colon_pos != std::string::npos) {
        host = remaining.substr(0, colon_pos);
        
        size_t port_end = (slash_pos != std::string::npos) ? slash_pos : remaining.length();
        std::string port_str = remaining.substr(colon_pos + 1, port_end - colon_pos - 1);
        
        try {
            port = std::stoi(port_str);
        } catch (...) {
            port = 6379; // Default
        }
    } else {
        if (slash_pos != std::string::npos) {
            host = remaining.substr(0, slash_pos);
        } else {
            host = remaining;
        }
    }
    
    // Parse database number
    if (slash_pos != std::string::npos && slash_pos + 1 < remaining.length()) {
        std::string db_str = remaining.substr(slash_pos + 1);
        try {
            db = std::stoi(db_str);
        } catch (...) {
            db = 0; // Default
        }
    }
}

bool RedisClient::testConnection() {
    if (!context_ || context_->err) {
        return false;
    }
    
    redisReply* reply = (redisReply*)redisCommand(context_.get(), "PING");
    if (reply == nullptr) {
        return false;
    }
    
    bool result = (reply->type == REDIS_REPLY_STATUS && 
                   strcmp(reply->str, "PONG") == 0);
    freeReplyObject(reply);
    return result;
}

bool RedisClient::pushJob(const std::string& simulation_id) {
    if (!context_ || context_->err) {
        return false;
    }
    
    redisReply* reply = (redisReply*)redisCommand(context_.get(), 
                                                   "LPUSH %s %s", 
                                                   queue_name_.c_str(), 
                                                   simulation_id.c_str());
    if (reply == nullptr) {
        return false;
    }
    
    bool result = (reply->type == REDIS_REPLY_INTEGER);
    freeReplyObject(reply);
    return result;
}

std::string RedisClient::popJob(int timeout_seconds) {
    if (!context_ || context_->err) {
        return "";
    }
    
    std::string result;
    
    if (timeout_seconds > 0) {
        // Blocking pop with timeout
        redisReply* reply = (redisReply*)redisCommand(context_.get(), 
                                                      "BRPOP %s %d", 
                                                      queue_name_.c_str(), 
                                                      timeout_seconds);
        if (reply == nullptr || reply->type != REDIS_REPLY_ARRAY || reply->elements < 2) {
            if (reply) freeReplyObject(reply);
            return "";
        }
        
        // Reply is an array: [queue_name, value]
        if (reply->element[1]->type == REDIS_REPLY_STRING) {
            result = std::string(reply->element[1]->str, reply->element[1]->len);
        }
        freeReplyObject(reply);
    } else {
        // Non-blocking pop
        result = tryPopJob();
    }
    
    return result;
}

std::string RedisClient::tryPopJob() {
    if (!context_ || context_->err) {
        return "";
    }
    
    redisReply* reply = (redisReply*)redisCommand(context_.get(), 
                                                   "RPOP %s", 
                                                   queue_name_.c_str());
    if (reply == nullptr || reply->type != REDIS_REPLY_STRING) {
        if (reply) freeReplyObject(reply);
        return "";
    }
    
    std::string result(reply->str, reply->len);
    freeReplyObject(reply);
    return result;
}

int RedisClient::getQueueLength() {
    if (!context_ || context_->err) {
        return -1;
    }
    
    redisReply* reply = (redisReply*)redisCommand(context_.get(), 
                                                   "LLEN %s", 
                                                   queue_name_.c_str());
    if (reply == nullptr || reply->type != REDIS_REPLY_INTEGER) {
        if (reply) freeReplyObject(reply);
        return -1;
    }
    
    int length = static_cast<int>(reply->integer);
    freeReplyObject(reply);
    return length;
}

