import { NextRequest, NextResponse } from 'next/server';
import { getSimulationMetrics, getSimulation } from '@/lib/simulations';

// GET /api/v1/simulations/[id]/metrics - Get simulation metrics
export async function GET(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> }
) {
  try {
    const { id } = await params;
    
    // Verify simulation exists
    const simulation = await getSimulation(id);
    if (!simulation) {
      return NextResponse.json(
        { error: 'Simulation not found' },
        { status: 404 }
      );
    }
    
    // Get metrics
    const metrics = await getSimulationMetrics(id);
    
    // Transform metrics to match the format expected by the frontend
    // Group by metric type for easier consumption
    const transformed = {
      euler_angles: metrics.map(m => m.euler_angles),
      ang_mom_body_frame: metrics.map(m => m.ang_mom_body_frame),
      a_control_torque: metrics.map(m => m.a_control_torque),
      a_command: metrics.map(m => m.a_command),
      state: metrics.map(m => m.state),
      d: metrics.map(m => [m.distance]), // Wrap in array to match expected format
      t: metrics.map(m => [m.timestamp]), // Wrap in array to match expected format
    };
    
    return NextResponse.json(transformed);
  } catch (error) {
    console.error('Error getting simulation metrics:', error);
    return NextResponse.json(
      { error: 'Failed to get simulation metrics' },
      { status: 500 }
    );
  }
}

