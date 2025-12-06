import { NextRequest, NextResponse } from 'next/server';
import { getSimulationMetrics, getSimulation } from '@/lib/simulations';

// GET /api/v1/simulations/[id]/metrics - Get simulation metrics
// Query parameters:
//   - since_step: Get metrics with step_index > since_step (optional)
//   - limit: Limit the number of results (optional)
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
    
    // Parse query parameters
    const searchParams = request.nextUrl.searchParams;
    const sinceStepParam = searchParams.get('since_step');
    const limitParam = searchParams.get('limit');
    
    const sinceStep = sinceStepParam ? parseInt(sinceStepParam, 10) : undefined;
    const limit = limitParam ? parseInt(limitParam, 10) : undefined;
    
    // Validate parameters
    if (sinceStepParam !== null && (isNaN(sinceStep!) || sinceStep! < 0)) {
      return NextResponse.json(
        { error: 'Invalid since_step parameter. Must be a non-negative integer.' },
        { status: 400 }
      );
    }
    
    if (limitParam !== null && (isNaN(limit!) || limit! <= 0)) {
      return NextResponse.json(
        { error: 'Invalid limit parameter. Must be a positive integer.' },
        { status: 400 }
      );
    }
    
    // Get metrics
    const metrics = await getSimulationMetrics(id, sinceStep, limit);
    
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
      t: [metrics.map(m => m.timestamp)], // Single series, wrap in array
    };
    
    // Get the last step_index to help client track progress
    const lastStep = metrics.length > 0 ? metrics[metrics.length - 1].step_index : (sinceStep ?? -1);
    
    // Check if there might be more data (if limit was applied and we got that many results)
    const hasMore = limit !== undefined && metrics.length === limit;
    
    return NextResponse.json({
      ...transformed,
      metadata: {
        count: metrics.length,
        lastStep,
        hasMore: hasMore,
      },
    });
  } catch (error) {
    console.error('Error getting simulation metrics:', error);
    return NextResponse.json(
      { error: 'Failed to get simulation metrics' },
      { status: 500 }
    );
  }
}

