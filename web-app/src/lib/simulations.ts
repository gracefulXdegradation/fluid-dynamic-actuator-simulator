import { prisma } from './db';
import { SimulationStatus as PrismaSimulationStatus, Prisma } from '@prisma/client';

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
  const simulation = await prisma.simulation.create({
    data: {
      status: 'scheduled',
      input_parameters: inputParams as unknown as Prisma.InputJsonValue,
    },
  });
  
  return simulation.id;
}

// Get a simulation by ID
export async function getSimulation(id: string): Promise<Simulation | null> {
  const simulation = await prisma.simulation.findUnique({
    where: { id },
  });
  
  if (!simulation) {
    return null;
  }
  
  return {
    id: simulation.id,
    status: simulation.status as SimulationStatus,
    created_at: simulation.created_at,
    started_at: simulation.started_at,
    completed_at: simulation.completed_at,
    input_parameters: simulation.input_parameters as unknown as SimulationInputParams,
    error_message: simulation.error_message,
  };
}

// List all simulations with optional status filter
export async function listSimulations(
  status?: SimulationStatus
): Promise<Simulation[]> {
  const simulations = await prisma.simulation.findMany({
    where: status ? { status: status as PrismaSimulationStatus } : undefined,
    orderBy: { created_at: 'desc' },
  });
  
  return simulations.map((simulation) => ({
    id: simulation.id,
    status: simulation.status as SimulationStatus,
    created_at: simulation.created_at,
    started_at: simulation.started_at,
    completed_at: simulation.completed_at,
    input_parameters: simulation.input_parameters as unknown as SimulationInputParams,
    error_message: simulation.error_message,
  }));
}

// Get metrics for a simulation
export async function getSimulationMetrics(
  simulationId: string,
  sinceStep?: number,
  limit?: number
): Promise<SimulationMetric[]> {
  const where: Prisma.SimulationMetricWhereInput = {
    simulation_id: simulationId,
    ...(sinceStep !== undefined && { step_index: { gt: sinceStep } }),
  };
  
  const metrics = await prisma.simulationMetric.findMany({
    where,
    orderBy: { step_index: 'asc' },
    take: limit, // limit is undefined if not provided, which means no limit
  });
  
  return metrics.map((metric) => ({
    id: metric.id,
    simulation_id: metric.simulation_id,
    step_index: metric.step_index,
    timestamp: Number(metric.timestamp), // Convert BigInt to number
    euler_angles: metric.euler_angles,
    ang_mom_body_frame: metric.ang_mom_body_frame,
    a_control_torque: metric.a_control_torque,
    a_command: metric.a_command,
    state: metric.state,
    distance: metric.distance,
  }));
}

// Delete a simulation (and its metrics via CASCADE)
export async function deleteSimulation(id: string): Promise<boolean> {
  try {
    await prisma.simulation.delete({
      where: { id },
    });
    return true;
  } catch (error) {
    // Prisma throws P2025 if record not found
    if (error && typeof error === 'object' && 'code' in error && error.code === 'P2025') {
      return false;
    }
    throw error;
  }
}

