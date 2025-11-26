import { NextRequest, NextResponse } from 'next/server';
import { createSimulation, listSimulations, SimulationInputParams, SimulationStatus } from '@/lib/simulations';
import { RedisClient } from '@/lib/redis';

// GET /api/v1/simulations - List all simulations
export async function GET(request: NextRequest) {
  try {
    const searchParams = request.nextUrl.searchParams;
    const status = searchParams.get('status') as SimulationStatus | undefined;
    
    const simulations = await listSimulations(status || undefined);
    
    return NextResponse.json({ simulations });
  } catch (error) {
    console.error('Error listing simulations:', error);
    return NextResponse.json(
      { error: 'Failed to list simulations' },
      { status: 500 }
    );
  }
}

// POST /api/v1/simulations - Create a new simulation
export async function POST(request: NextRequest) {
  try {
    const body = await request.json();
    
    // Convert datetime-local format (YYYY-MM-DDTHH:mm) to required format (YYYY-MM-DD HH:mm:ss)
    const formatDateTime = (dateTimeStr: string): string => {
      if (!dateTimeStr) return dateTimeStr;
      // Replace 'T' with ' ' and add ':00' for seconds
      return dateTimeStr.replace('T', ' ') + ':00';
    };
    
    // Validate input parameters
    const inputParams: SimulationInputParams = {
      tle_line1: body.tle_line1,
      tle_line2: body.tle_line2,
      start_date_time: formatDateTime(body.start_date_time),
      end_date_time: formatDateTime(body.end_date_time),
      control_time_step: body.control_time_step,
      ground_station_lla: {
        lat: body.ground_station_lla?.lat,
        lon: body.ground_station_lla?.lon,
        alt: body.ground_station_lla?.alt,
      },
      ground_station_elevation: body.ground_station_elevation,
    };
    
    // Basic validation
    if (!inputParams.tle_line1 || !inputParams.tle_line2) {
      return NextResponse.json(
        { error: 'TLE lines are required' },
        { status: 400 }
      );
    }
    
    if (!inputParams.start_date_time || !inputParams.end_date_time) {
      return NextResponse.json(
        { error: 'Start and end date times are required' },
        { status: 400 }
      );
    }
    
    if (!inputParams.control_time_step || inputParams.control_time_step <= 0) {
      return NextResponse.json(
        { error: 'Control time step must be a positive number' },
        { status: 400 }
      );
    }
    
    // Create simulation in database
    const simulationId = await createSimulation(inputParams);
    
    // Push job to Redis queue
    try {
      const redis = new RedisClient();
      await redis.pushJob(simulationId);
    } catch (redisError) {
      console.error('Error pushing job to Redis:', redisError);
      // Continue even if Redis fails - the simulation is created and can be processed later
    }
    
    return NextResponse.json({ 
      id: simulationId,
      message: 'Simulation created and queued'
    }, { status: 201 });
  } catch (error) {
    console.error('Error creating simulation:', error);
    return NextResponse.json(
      { error: 'Failed to create simulation' },
      { status: 500 }
    );
  }
}

