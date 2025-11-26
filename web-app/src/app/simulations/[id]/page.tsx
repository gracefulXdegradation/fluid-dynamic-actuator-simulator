"use client"
import { useState, useEffect } from 'react';
import { useParams } from 'next/navigation'
import LineGraph from '@/components/LineGraph';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Skeleton } from '@/components/ui/skeleton';
import { Alert, AlertDescription } from '@/components/ui/alert';

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

const SimulationPage = () => {
  const [data, setData] = useState<Data | null>(null);
  const [loading, setLoading] = useState(true);
  const params = useParams<{ id: string }>()
  const { id } = params;

  useEffect(() => {
    if (id) {
      const fetchData = async () => {
        try {
          const response = await fetch(`/api/v1/simulations/${id}/metrics`);
          if (!response.ok) {
            throw new Error('Failed to fetch simulation metrics');
          }
          const result = await response.json();
          setData(result);
        } catch (error) {
          console.error('Error fetching data:', error);
        } finally {
          setLoading(false);
        }
      };
      fetchData();
    }
  }, [id]);

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

  if (!data) {
    return (
      <div className="page-container">
        <div className="w-full max-w-7xl mx-auto">
          <h1 className="text-2xl font-bold mb-6">Simulation #{id}</h1>
          <Alert variant="destructive">
            <AlertDescription>No data available for this simulation.</AlertDescription>
          </Alert>
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
        <h1 className="text-2xl font-bold mb-6">Simulation #{id}</h1>
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
