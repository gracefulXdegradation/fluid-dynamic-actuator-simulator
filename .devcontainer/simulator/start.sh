#!/bin/bash
set -e

echo "========================================="
echo "Starting simulator container..."
echo "========================================="

# Wait for PostgreSQL to be ready
echo "Waiting for PostgreSQL..."
until pg_isready -h postgres -U fds_user -d fds_db > /dev/null 2>&1; do
  echo "PostgreSQL is unavailable - sleeping"
  sleep 1
done
echo "✓ PostgreSQL is ready!"

# Wait for Redis to be ready
echo "Waiting for Redis..."
until redis-cli -h redis ping > /dev/null 2>&1; do
  echo "Redis is unavailable - sleeping"
  sleep 1
done
echo "✓ Redis is ready!"

# Change to simulator directory
cd /workspace/simulator

# Build the project if build directory doesn't exist or job_processor doesn't exist
if [ ! -d "build" ] || [ ! -f "build/bin/job_processor" ]; then
  echo "Building C++ targets..."
  cmake -S . -B build
  cmake --build build
  echo "✓ Build complete!"
else
  echo "Build directory exists, checking if rebuild is needed..."
  # Optionally, you could add logic here to check if source files changed
  # For now, we'll just use the existing build
fi

# Verify job processor exists
if [ ! -f "build/bin/job_processor" ]; then
  echo "ERROR: job_processor executable not found after build!"
  exit 1
fi

# Start the job processor
echo "========================================="
echo "Starting job processor..."
echo "========================================="
exec /workspace/simulator/build/bin/job_processor

