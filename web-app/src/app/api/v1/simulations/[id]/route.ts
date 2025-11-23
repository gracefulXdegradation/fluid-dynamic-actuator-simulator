import { NextRequest, NextResponse } from 'next/server';
import { getSimulation, getSimulationMetrics, deleteSimulation } from '@/lib/simulations';

// GET /api/v1/simulations/[id] - Get simulation details
export async function GET(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> }
) {
  try {
    const { id } = await params;
    const simulation = await getSimulation(id);
    
    if (!simulation) {
      return NextResponse.json(
        { error: 'Simulation not found' },
        { status: 404 }
      );
    }
    
    return NextResponse.json({ simulation });
  } catch (error) {
    console.error('Error getting simulation:', error);
    return NextResponse.json(
      { error: 'Failed to get simulation' },
      { status: 500 }
    );
  }
}

// DELETE /api/v1/simulations/[id] - Delete simulation
export async function DELETE(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> }
) {
  try {
    const { id } = await params;
    const deleted = await deleteSimulation(id);
    
    if (!deleted) {
      return NextResponse.json(
        { error: 'Simulation not found' },
        { status: 404 }
      );
    }
    
    return NextResponse.json({ message: 'Simulation deleted' });
  } catch (error) {
    console.error('Error deleting simulation:', error);
    return NextResponse.json(
      { error: 'Failed to delete simulation' },
      { status: 500 }
    );
  }
}

