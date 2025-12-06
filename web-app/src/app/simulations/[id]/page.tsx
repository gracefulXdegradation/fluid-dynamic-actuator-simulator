"use client"
import { useState, useEffect, useRef } from 'react';
import { useParams } from 'next/navigation'
import LineGraph from '@/components/LineGraph';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Skeleton } from '@/components/ui/skeleton';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Badge } from '@/components/ui/badge';

const rad2deg = (rad: number) => rad * 180 / Math.PI;

interface Data {
  euler_angles: number[][];
  ang_mom_body_frame: number[][];
  a_control_torque: number[][];
  a_command: number[][];
  state: number[][];
  d: number[][];
  t: number[][];
}

type SimulationStatus = 'scheduled' | 'running' | 'completed' | 'failed';
type ConnectionStatus = 'disconnected' | 'connecting' | 'connected';

const SimulationPage = () => {
  const [data, setData] = useState<Data | null>(null);
  const [loading, setLoading] = useState(true);
  const [simulationStatus, setSimulationStatus] = useState<SimulationStatus | null>(null);
  const [connectionStatus, setConnectionStatus] = useState<ConnectionStatus>('disconnected');
  const eventSourceRef = useRef<EventSource | null>(null);
  const params = useParams<{ id: string }>()
  const { id } = params;

  // Helper function to merge new metrics with existing data
  const mergeMetrics = (existing: Data | null, newMetrics: Data): Data => {
    if (!existing) {
      return newMetrics;
    }

    // Merge arrays by concatenating
    return {
      euler_angles: existing.euler_angles.map((series, idx) => [...series, ...newMetrics.euler_angles[idx]]),
      ang_mom_body_frame: existing.ang_mom_body_frame.map((series, idx) => [...series, ...newMetrics.ang_mom_body_frame[idx]]),
      a_control_torque: existing.a_control_torque.map((series, idx) => [...series, ...newMetrics.a_control_torque[idx]]),
      a_command: existing.a_command.map((series, idx) => [...series, ...newMetrics.a_command[idx]]),
      state: existing.state.map((series, idx) => [...series, ...newMetrics.state[idx]]),
      d: existing.d.map((series, idx) => [...series, ...newMetrics.d[idx]]),
      t: existing.t.map((series, idx) => [...series, ...newMetrics.t[idx]]),
    };
  };

  useEffect(() => {
    if (!id) return;

    let isMounted = true;

    // Fetch simulation status first
    const fetchSimulationStatus = async () => {
      try {
        const response = await fetch(`/api/v1/simulations/${id}`);
        if (!response.ok) {
          throw new Error('Failed to fetch simulation');
        }
        const result = await response.json();
        const status = result.simulation?.status as SimulationStatus;
        
        if (isMounted) {
          setSimulationStatus(status);
          
          // If simulation is completed or failed, fetch all metrics once
          if (status === 'completed' || status === 'failed') {
            try {
              const metricsResponse = await fetch(`/api/v1/simulations/${id}/metrics`);
              if (metricsResponse.ok) {
                const metricsResult = await metricsResponse.json();
                // Extract only the data fields, ignoring metadata if present
                const metricsData: Data = {
                  euler_angles: metricsResult.euler_angles,
                  ang_mom_body_frame: metricsResult.ang_mom_body_frame,
                  a_control_torque: metricsResult.a_control_torque,
                  a_command: metricsResult.a_command,
                  state: metricsResult.state,
                  d: metricsResult.d,
                  t: metricsResult.t,
                };
                setData(metricsData);
              }
            } catch (error) {
              console.error('Error fetching completed simulation metrics:', error);
            } finally {
              if (isMounted) {
                setLoading(false);
              }
            }
          } else if (status === 'running' || status === 'scheduled') {
            // For running/scheduled simulations, set up SSE connection
            setLoading(false); // Show the page even if no data yet
            setupSSEConnection();
          } else {
            setLoading(false);
          }
        }
      } catch (error) {
        console.error('Error fetching simulation status:', error);
        if (isMounted) {
          setLoading(false);
        }
      }
    };

    // Set up SSE connection for real-time updates
    const setupSSEConnection = () => {
      if (eventSourceRef.current) {
        eventSourceRef.current.close();
      }

      setConnectionStatus('connecting');
      const eventSource = new EventSource(`/api/v1/simulations/${id}/metrics/stream`);
      eventSourceRef.current = eventSource;

      eventSource.onopen = () => {
        if (isMounted) {
          setConnectionStatus('connected');
        }
      };

      eventSource.onmessage = (event) => {
        try {
          const message = JSON.parse(event.data);
          
          if (message.type === 'connected') {
            if (isMounted) {
              setConnectionStatus('connected');
            }
          } else if (message.type === 'metrics') {
            if (isMounted && message.data) {
              setData((prevData) => mergeMetrics(prevData, message.data));
            }
          } else if (message.type === 'status') {
            if (isMounted) {
              setSimulationStatus(message.status);
              if (message.status === 'completed' || message.status === 'failed') {
                eventSource.close();
                setConnectionStatus('disconnected');
              }
            }
          } else if (message.type === 'error') {
            console.error('SSE error:', message.message);
          } else if (message.type === 'keepalive') {
            // Keepalive - no action needed
          }
        } catch (error) {
          console.error('Error parsing SSE message:', error);
        }
      };

      eventSource.onerror = (error) => {
        console.error('SSE connection error:', error);
        if (isMounted) {
          setConnectionStatus('disconnected');
          // Retry connection after 3 seconds
          setTimeout(() => {
            if (isMounted && (simulationStatus === 'running' || simulationStatus === 'scheduled')) {
              setupSSEConnection();
            }
          }, 3000);
        }
      };
    };

    fetchSimulationStatus();

    // Cleanup on unmount
    return () => {
      isMounted = false;
      if (eventSourceRef.current) {
        eventSourceRef.current.close();
        eventSourceRef.current = null;
      }
    };
  }, [id, simulationStatus]);

  if (loading) {
    return (
      <div className="page-container">
        <div className="w-full max-w-7xl mx-auto">
          <Skeleton className="h-8 w-32 mb-6" />
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            {[1, 2, 3, 4, 5, 6, 7].map((i) => (
              <Card key={i}>
                <CardHeader>
                  <Skeleton className="h-6 w-48" />
                </CardHeader>
                <CardContent>
                  <Skeleton className="h-64 w-full" />
                </CardContent>
              </Card>
            ))}
          </div>
        </div>
      </div>
    );
  }

  // Get connection status badge
  const getConnectionBadge = () => {
    if (simulationStatus === 'running' || simulationStatus === 'scheduled') {
      switch (connectionStatus) {
        case 'connecting':
          return <Badge variant="outline" className="ml-2">Connecting...</Badge>;
        case 'connected':
          return <Badge variant="default" className="ml-2 bg-green-500">Live</Badge>;
        case 'disconnected':
          return <Badge variant="outline" className="ml-2">Reconnecting...</Badge>;
      }
    }
    return null;
  };

  if (!data) {
    const isRunning = simulationStatus === 'running' || simulationStatus === 'scheduled';
    return (
      <div className="page-container">
        <div className="w-full max-w-7xl mx-auto">
          <div className="flex items-center mb-6">
            <h1 className="text-2xl font-bold">Simulation #{id}</h1>
            {getConnectionBadge()}
          </div>
          {isRunning ? (
            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
              {[1, 2, 3, 4, 5, 6, 7].map((i) => (
                <Card key={i}>
                  <CardHeader>
                    <Skeleton className="h-6 w-48" />
                  </CardHeader>
                  <CardContent>
                    <Skeleton className="h-64 w-full" />
                  </CardContent>
                </Card>
              ))}
            </div>
          ) : (
            <Alert variant="destructive">
              <AlertDescription>No data available for this simulation.</AlertDescription>
            </Alert>
          )}
        </div>
      </div>
    );
  }

  // Prepare data for plotting
  const { t, a_control_torque, a_command, ang_mom_body_frame, euler_angles, state, d } = data;
  
  const angularRate = state.slice(4,7).map((data: number[]) => data.map(rad2deg));
  const angularRateAbs = angularRate[0].map((_, i) => Math.pow(angularRate[0][i], 2) + Math.pow(angularRate[1][i], 2) + Math.pow(angularRate[2][i],2) )
  const angularMomentum = [7, 9, 11 ,13].map(i => state[i].map(v => v * 1e6));
  const angularMomentumBodyFrame = ang_mom_body_frame.map(dim => dim.map(v => v * 1e6));
  const ts = t[0];

  return (
    <div className="page-container">
      <div className="w-full max-w-7xl mx-auto">
        <div className="flex items-center mb-6">
          <h1 className="text-2xl font-bold">Simulation #{id}</h1>
          {getConnectionBadge()}
        </div>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <Card>
            <CardHeader>
              <CardTitle className="text-lg">Required control torque</CardTitle>
            </CardHeader>
            <CardContent>
              <LineGraph
                timestamps={ts}
                values={a_control_torque}
                graphNames={["&tau;<sub>1</sub>", "&tau;<sub>2</sub>", "&tau;<sub>3</sub>", "&tau;<sub>4</sub>"]}
                labelX="Time"
                labelY="Torque [mNm]"
              />
            </CardContent>
          </Card>
          <Card>
            <CardHeader>
              <CardTitle className="text-lg">Body angular rate w.r.t. body frame</CardTitle>
            </CardHeader>
            <CardContent>
              <LineGraph
                timestamps={ts}
                values={[...angularRate, angularRateAbs]}
                graphNames={["<sub>b</sub>&omega;<sub>bx</sub>", "<sub>b</sub>&omega;<sub>by</sub>", "<sub>b</sub>&omega;<sub>bz</sub>", "|<sub>b</sub>&omega;<sub>b</sub>|"]}
                labelX="Time"
                labelY=""
              />
            </CardContent>
          </Card>
          <Card>
            <CardHeader>
              <CardTitle className="text-lg">Actuator angular momentum in actuator frame</CardTitle>
            </CardHeader>
            <CardContent>
              <LineGraph
                timestamps={ts}
                values={angularMomentum}
                graphNames={["h<sub>1</sub>", "h<sub>2</sub>", "h<sub>3</sub>", "h<sub>4</sub>"]}
                labelX="Time"
                labelY="Angular momentum [&mu;Nms]"
              />
            </CardContent>
          </Card>
          <Card>
            <CardHeader>
              <CardTitle className="text-lg">Actuator angular momentum in body frame</CardTitle>
            </CardHeader>
            <CardContent>
              <LineGraph
                timestamps={ts}
                values={angularMomentumBodyFrame}
                graphNames={["h<sub>x</sub>", "h<sub>y</sub>", "h<sub>z</sub>"]}
                labelX="Time"
                labelY="Angular momentum [&mu;Nms]"
              />
            </CardContent>
          </Card>
          <Card>
            <CardHeader>
              <CardTitle className="text-lg">Distance</CardTitle>
            </CardHeader>
            <CardContent>
              <LineGraph
                timestamps={ts}
                values={d}
                graphNames={["Distance"]}
                labelX="Time"
                labelY="Distance to the ground station [km]"
              />
            </CardContent>
          </Card>
          <Card>
            <CardHeader>
              <CardTitle className="text-lg">Attitude error angle</CardTitle>
            </CardHeader>
            <CardContent>
              <LineGraph
                timestamps={ts}
                values={[euler_angles[1].map(rad => Math.min(Math.max(rad2deg(rad), 0), 0.2))]}
                graphNames={["y"]}
                labelX="Time"
                labelY="Error angle [&deg;]"
              />
            </CardContent>
          </Card>
          <Card>
            <CardHeader>
              <CardTitle className="text-lg">Actuator commands</CardTitle>
            </CardHeader>
            <CardContent>
              <LineGraph
                timestamps={ts}
                values={a_command}
                graphNames={["&mu;<sub>1</sub>", "&mu;<sub>2</sub>", "&mu;<sub>3</sub>", "&mu;<sub>4</sub>"]}
                labelX="Time"
                labelY=""
              />
            </CardContent>
          </Card>
        </div>
      </div>
    </div>
  );
};

export default SimulationPage;
