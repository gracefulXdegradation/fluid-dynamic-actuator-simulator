import { NextRequest } from 'next/server';
import { getSimulationMetrics, getSimulation } from '@/lib/simulations';

// GET /api/v1/simulations/[id]/metrics/stream - Stream simulation metrics via Server-Sent Events
export async function GET(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> }
) {
  const { id } = await params;
  
  // Verify simulation exists
  const simulation = await getSimulation(id);
  if (!simulation) {
    return new Response(
      JSON.stringify({ error: 'Simulation not found' }),
      {
        status: 404,
        headers: { 'Content-Type': 'application/json' },
      }
    );
  }

  // Create a readable stream for SSE
  const stream = new ReadableStream({
    async start(controller) {
      const encoder = new TextEncoder();
      let lastStep = -1;
      let isActive = true;
      
      // Helper function to send SSE event
      const sendEvent = (type: string, data: any) => {
        if (!isActive) return;
        try {
          const message = `data: ${JSON.stringify({ type, ...data })}\n\n`;
          controller.enqueue(encoder.encode(message));
        } catch (error) {
          console.error('Error sending SSE event:', error);
        }
      };

      // Helper function to transform metrics to frontend format
      const transpose = (arr: number[][]): number[][] => {
        if (arr.length === 0) return [];
        const numSeries = arr[0].length;
        const result: number[][] = [];
        for (let i = 0; i < numSeries; i++) {
          result.push(arr.map(row => row[i]));
        }
        return result;
      };

      // Poll for new metrics
      const pollInterval = setInterval(async () => {
        if (!isActive) {
          clearInterval(pollInterval);
          return;
        }

        try {
          // Check simulation status
          const currentSimulation = await getSimulation(id);
          if (!currentSimulation) {
            sendEvent('error', { message: 'Simulation not found' });
            clearInterval(pollInterval);
            controller.close();
            return;
          }

          // If simulation is completed or failed, send final status and close
          if (currentSimulation.status === 'completed' || currentSimulation.status === 'failed') {
            // Get any remaining metrics
            const remainingMetrics = await getSimulationMetrics(id, lastStep);
            
            if (remainingMetrics.length > 0) {
              const transformed = {
                euler_angles: transpose(remainingMetrics.map(m => m.euler_angles)),
                ang_mom_body_frame: transpose(remainingMetrics.map(m => m.ang_mom_body_frame)),
                a_control_torque: transpose(remainingMetrics.map(m => m.a_control_torque)),
                a_command: transpose(remainingMetrics.map(m => m.a_command)),
                state: transpose(remainingMetrics.map(m => m.state)),
                d: [remainingMetrics.map(m => m.distance)],
                t: [remainingMetrics.map(m => m.timestamp)],
              };
              
              sendEvent('metrics', {
                data: transformed,
                lastStep: remainingMetrics[remainingMetrics.length - 1].step_index,
              });
              lastStep = remainingMetrics[remainingMetrics.length - 1].step_index;
            }
            
            sendEvent('status', { status: currentSimulation.status });
            clearInterval(pollInterval);
            controller.close();
            return;
          }

          // Get new metrics since last step
          const newMetrics = await getSimulationMetrics(id, lastStep);
          
          if (newMetrics.length > 0) {
            const transformed = {
              euler_angles: transpose(newMetrics.map(m => m.euler_angles)),
              ang_mom_body_frame: transpose(newMetrics.map(m => m.ang_mom_body_frame)),
              a_control_torque: transpose(newMetrics.map(m => m.a_control_torque)),
              a_command: transpose(newMetrics.map(m => m.a_command)),
              state: transpose(newMetrics.map(m => m.state)),
              d: [newMetrics.map(m => m.distance)],
              t: [newMetrics.map(m => m.timestamp)],
            };
            
            sendEvent('metrics', {
              data: transformed,
              lastStep: newMetrics[newMetrics.length - 1].step_index,
            });
            
            lastStep = newMetrics[newMetrics.length - 1].step_index;
          }
        } catch (error) {
          console.error('Error polling metrics:', error);
          sendEvent('error', { message: 'Failed to fetch metrics' });
        }
      }, 2000); // Poll every 2 seconds

      // Send initial connection message
      sendEvent('connected', { simulationId: id });

      // Handle client disconnect
      request.signal.addEventListener('abort', () => {
        isActive = false;
        clearInterval(pollInterval);
        controller.close();
      });

      // Send keepalive every 30 seconds to maintain connection
      const keepaliveInterval = setInterval(() => {
        if (!isActive) {
          clearInterval(keepaliveInterval);
          return;
        }
        sendEvent('keepalive', {});
      }, 30000);
    },
  });

  // Return SSE response
  return new Response(stream, {
    headers: {
      'Content-Type': 'text/event-stream',
      'Cache-Control': 'no-cache, no-transform',
      'Connection': 'keep-alive',
      'X-Accel-Buffering': 'no', // Disable buffering in nginx
    },
  });
}

