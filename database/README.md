# Database Schema

This directory contains database initialization scripts that run automatically when the PostgreSQL container is first created.

## Schema Overview

### `simulations` table
Stores simulation jobs with their input parameters and execution status.

- `id`: UUID primary key
- `status`: Enum (scheduled, running, completed, failed)
- `created_at`: Timestamp when simulation was created
- `started_at`: Timestamp when simulation started running
- `completed_at`: Timestamp when simulation completed
- `input_parameters`: JSONB object containing:
  - `tle_line1`: First line of TLE data
  - `tle_line2`: Second line of TLE data
  - `start_date_time`: Simulation start time (ISO 8601 format)
  - `end_date_time`: Simulation end time (ISO 8601 format)
  - `control_time_step`: Time step in milliseconds
  - `ground_station_lla`: Object with `lat`, `lon`, `alt`
  - `ground_station_elevation`: Elevation angle in degrees
- `error_message`: Error message if simulation failed

### `simulation_metrics` table
Stores metrics for each time step of a simulation.

- `id`: UUID primary key
- `simulation_id`: Foreign key to simulations table
- `step_index`: Zero-based index of the time step
- `timestamp`: Unix timestamp in milliseconds
- `euler_angles`: Array of 3 Euler angles (radians)
- `ang_mom_body_frame`: Array of 3 angular momentum components in body frame
- `a_control_torque`: Array of 4 actuator control torques
- `a_command`: Array of 4 actuator commands
- `state`: Array of 15 state variables (quaternion, angular rate, actuator states)
- `distance`: Distance to ground station (meters)

## Connection Details

- Host: `postgres` (from within Docker network) or `localhost` (from host)
- Port: `5432`
- Database: `fds_db`
- User: `fds_user`
- Password: `fds_password`

## Environment Variables

Both `simulator` and `web-app` containers have access to:
- `DATABASE_URL`: `postgresql://fds_user:fds_password@postgres:5432/fds_db`
- `REDIS_URL`: `redis://redis:6379`

