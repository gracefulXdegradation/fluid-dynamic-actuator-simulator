-- Create enum type for simulation status
CREATE TYPE simulation_status AS ENUM ('scheduled', 'running', 'completed', 'failed');

-- Create simulations table
CREATE TABLE simulations (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    status simulation_status NOT NULL DEFAULT 'scheduled',
    created_at TIMESTAMP WITH TIME ZONE NOT NULL DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP WITH TIME ZONE,
    completed_at TIMESTAMP WITH TIME ZONE,
    input_parameters JSONB NOT NULL,
    error_message TEXT
);

-- Create simulation_metrics table
CREATE TABLE simulation_metrics (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    simulation_id UUID NOT NULL REFERENCES simulations(id) ON DELETE CASCADE,
    step_index INTEGER NOT NULL,
    timestamp BIGINT NOT NULL, -- milliseconds since epoch
    euler_angles DOUBLE PRECISION[3] NOT NULL,
    ang_mom_body_frame DOUBLE PRECISION[3] NOT NULL,
    a_control_torque DOUBLE PRECISION[4] NOT NULL,
    a_command DOUBLE PRECISION[4] NOT NULL,
    state DOUBLE PRECISION[15] NOT NULL,
    distance DOUBLE PRECISION NOT NULL,
    UNIQUE(simulation_id, step_index)
);

-- Create indexes for better query performance
CREATE INDEX idx_simulations_status ON simulations(status);
CREATE INDEX idx_simulations_created_at ON simulations(created_at DESC);
CREATE INDEX idx_simulation_metrics_simulation_id ON simulation_metrics(simulation_id);
CREATE INDEX idx_simulation_metrics_step_index ON simulation_metrics(simulation_id, step_index);

-- Add comments for documentation
COMMENT ON TABLE simulations IS 'Stores simulation jobs with their input parameters and status';
COMMENT ON TABLE simulation_metrics IS 'Stores metrics for each step of a simulation';
COMMENT ON COLUMN simulations.input_parameters IS 'JSON object containing: TLE lines (line1, line2), start_date_time, end_date_time, control_time_step, ground_station_lla, ground_station_elevation';
COMMENT ON COLUMN simulation_metrics.timestamp IS 'Unix timestamp in milliseconds since epoch';

