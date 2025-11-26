# Job Processor

The job processor is a daemon that polls Redis for simulation jobs and executes them automatically.

## Building

The job processor is built automatically when you build the simulator:

```bash
cd simulator
cmake -S . -B build
cmake --build build
```

This will create two executables:
- `build/bin/fds` - The main simulator (can be run directly with a simulation_id)
- `build/bin/job_processor` - The job processor daemon

## Running the Job Processor

### Prerequisites

The job processor requires the following environment variables:
- `DATABASE_URL` - PostgreSQL connection string (e.g., `postgresql://fds_user:fds_password@postgres:5432/fds_db`)
- `REDIS_URL` - Redis connection string (e.g., `redis://redis:6379`)

These are automatically set in the Docker container.

### Usage

```bash
# Run with default poll interval (5 seconds)
./build/bin/job_processor

# Run with custom poll interval (e.g., 10 seconds)
./build/bin/job_processor 10
```

The job processor will:
1. Connect to Redis and Database
2. Poll Redis for simulation jobs every N seconds (default: 5)
3. When a job is found, execute the simulator with that simulation_id
4. Continue polling until stopped (Ctrl+C or SIGTERM)

## How It Works

1. **Job Queue**: The web app pushes simulation IDs to a Redis list called `simulation_jobs`
2. **Polling**: The job processor polls Redis using `RPOP` (non-blocking)
3. **Execution**: When a job is found, it forks a process and executes `./build/bin/fds <simulation_id>`
4. **Status Updates**: The simulator updates the database status (running → completed/failed)

## Adding Jobs to the Queue

Jobs can be added to Redis from the web app or manually:

```bash
# Using redis-cli
docker exec -it fds-redis-1 redis-cli LPUSH simulation_jobs "your-simulation-id-here"

# Or using the RedisClient API from C++
RedisClient redis("redis://redis:6379");
redis.pushJob("your-simulation-id-here");
```

## Monitoring

You can check the queue length:

```bash
docker exec -it fds-redis-1 redis-cli LLEN simulation_jobs
```

You can view pending jobs (without removing them):

```bash
docker exec -it fds-redis-1 redis-cli LRANGE simulation_jobs 0 -1
```

## Stopping the Job Processor

- Press `Ctrl+C` to gracefully stop
- Or send `SIGTERM`: `kill <pid>`

The processor will finish processing the current job (if any) before stopping.

## Running Multiple Processors

You can run multiple job processor instances for parallel processing:

```bash
# Terminal 1
./build/bin/job_processor

# Terminal 2
./build/bin/job_processor

# Terminal 3
./build/bin/job_processor
```

Each processor will compete for jobs from the same Redis queue, allowing for parallel simulation execution.

