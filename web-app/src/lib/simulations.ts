import { query } from './db';
import { randomUUID } from 'crypto';

export type SimulationStatus = 'scheduled' | 'running' | 'completed' | 'failed';

export interface SimulationInputParams {
  tle_line1: string;
  tle_line2: string;
  start_date_time: string; // ISO 8601 format
  end_date_time: string; // ISO 8601 format
  control_time_step: number; // milliseconds
  ground_station_lla: {
    lat: number;
    lon: number;
    alt: number;
  };
  ground_station_elevation: number; // degrees
}

export interface Simulation {
  id: string;
  status: SimulationStatus;
  created_at: Date;
  started_at: Date | null;
  completed_at: Date | null;
  input_parameters: SimulationInputParams;
  error_message: string | null;
}

export interface SimulationMetric {
  id: string;
  simulation_id: string;
  step_index: number;
  timestamp: number; // milliseconds since epoch
  euler_angles: number[]; // [3]
  ang_mom_body_frame: number[]; // [3]
  a_control_torque: number[]; // [4]
  a_command: number[]; // [4]
  state: number[]; // [15]
  distance: number;
}

// Create a new simulation
export async function createSimulation(
  inputParams: SimulationInputParams
): Promise<string> {
  const simulationId = randomUUID();
  
  await query(
    `INSERT INTO simulations (id, status, input_parameters)
     VALUES ($1, $2, $3)`,
    [simulationId, 'scheduled', JSON.stringify(inputParams)]
  );
  
  return simulationId;
}

// Get a simulation by ID
export async function getSimulation(id: string): Promise<Simulation | null> {
  const result = await query<Simulation>(
    `SELECT 
      id,
      status,
      created_at,
      started_at,
      completed_at,
      input_parameters,
      error_message
     FROM simulations
     WHERE id = $1`,
    [id]
  );
  
  if (result.rows.length === 0) {
    return null;
  }
  
  const row = result.rows[0];
  return {
    ...row,
    input_parameters: row.input_parameters as SimulationInputParams,
  };
}

// List all simulations with optional status filter
export async function listSimulations(
  status?: SimulationStatus
): Promise<Simulation[]> {
  let sql = `SELECT 
    id,
    status,
    created_at,
    started_at,
    completed_at,
    input_parameters,
    error_message
   FROM simulations`;
  
  const params: any[] = [];
  
  if (status) {
    sql += ' WHERE status = $1';
    params.push(status);
  }
  
  sql += ' ORDER BY created_at DESC';
  
  const result = await query<Simulation>(sql, params);
  
  return result.rows.map((row) => ({
    ...row,
    input_parameters: row.input_parameters as SimulationInputParams,
  }));
}

// Get metrics for a simulation
export async function getSimulationMetrics(
  simulationId: string
): Promise<SimulationMetric[]> {
  const result = await query<SimulationMetric>(
    `SELECT 
      id,
      simulation_id,
      step_index,
      timestamp,
      euler_angles,
      ang_mom_body_frame,
      a_control_torque,
      a_command,
      state,
      distance
     FROM simulation_metrics
     WHERE simulation_id = $1
     ORDER BY step_index ASC`,
    [simulationId]
  );
  
  return result.rows;
}

// Delete a simulation (and its metrics via CASCADE)
export async function deleteSimulation(id: string): Promise<boolean> {
  const result = await query(
    'DELETE FROM simulations WHERE id = $1',
    [id]
  );
  
  return result.rowCount !== null && result.rowCount > 0;
}

