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
    
    // Helper function to transpose array of arrays
    // Converts from [[a0, b0, c0], [a1, b1, c1], ...] to [[a0, a1, ...], [b0, b1, ...], [c0, c1, ...]]
    const transpose = (arr: number[][]): number[][] => {
      if (arr.length === 0) return [];
      const numSeries = arr[0].length;
      const result: number[][] = [];
      for (let i = 0; i < numSeries; i++) {
        result.push(arr.map(row => row[i]));
      }
      return result;
    };
    
    // Transform metrics to match the format expected by the frontend
    // Convert from row-based (time steps) to column-based (series over time)
    const transformed = {
      euler_angles: transpose(metrics.map(m => m.euler_angles)),
      ang_mom_body_frame: transpose(metrics.map(m => m.ang_mom_body_frame)),
      a_control_torque: transpose(metrics.map(m => m.a_control_torque)),
      a_command: transpose(metrics.map(m => m.a_command)),
      state: transpose(metrics.map(m => m.state)),
      d: [metrics.map(m => m.distance)], // Single series, wrap in array
      t: [metrics.map(m => parseInt(m.timestamp, 10))], // Single series, wrap in array
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

