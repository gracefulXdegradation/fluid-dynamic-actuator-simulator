"use client"
import { useState, useEffect } from 'react';
import { useParams } from 'next/navigation'
import LineGraph from '@/components/LineGraph';

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
          const response = await fetch(`/api/v1/entries/${id}`);
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
      <div className="loading-container">
        <p>Loading data...</p>
      </div>
    );
  }

  if (!data) {
    return <p>No data available for this simulation.</p>;
  }

  // Prepare data for plotting
  const { t, a_control_torque, a_command, ang_mom_body_frame, euler_angles, state, d } = data;

  
  const angularRate = state.slice(4,7).map((data: number[]) => data.map(rad2deg));
  const angularRateAbs = angularRate[0].map((_, i) => Math.pow(angularRate[0][i], 2) + Math.pow(angularRate[1][i], 2) + Math.pow(angularRate[2][i],2) )
  const angularMomentum = [7, 9, 11 ,13].map(i => state[i].map(v => v * 1e6));
  const angularMomentumBodyFrame = ang_mom_body_frame.map(dim => dim.map(v => v * 1e6));
  const ts = t[0];

  return (
    <>
      <h1>Simulation {id}</h1>
      <div className='container'>
        <div>
          <h2>Required control torque</h2>
          <LineGraph
            timestamps={ts}
            values={a_control_torque}
            graphNames={["Tau_1", "Tau_2", "Tau_3", "Tau_4"]}
          />
        </div>
        <div>
          <h2>Body angular rate w.r.t. body frame</h2>
          <LineGraph
            timestamps={ts}
            values={[...angularRate, angularRateAbs]}
            graphNames={["Omega x", "Omega y", "Omega z", "|Omega|"]}
          />
        </div>
        <div>
          <h2>Actuator angular momentum in actuator frame</h2>
          <LineGraph
            timestamps={ts}
            values={angularMomentum}
            graphNames={["h1", "h2", "h3", "h4"]}
          />
        </div>
        <div>
          <h2>Actuator angular momentum in body frame</h2>
          <LineGraph
            timestamps={ts}
            values={angularMomentumBodyFrame}
            graphNames={["h_x", "h_y", "h_z"]}
          />
        </div>
        <div>
          <h2>Distance</h2>
          <LineGraph
            timestamps={ts}
            values={d}
            graphNames={["Distance"]}
          />
        </div>
        <div>
          <h2>Attitude error angle</h2>
          <LineGraph
            timestamps={ts}
            values={[euler_angles[1].map(rad => Math.min(Math.max(rad2deg(rad), 0), 0.2))]}
            graphNames={["y"]}
          />
        </div>
        <div>
          <h2>Actuator commands</h2>
          <LineGraph
            timestamps={ts}
            values={a_command}
            graphNames={["u1", "u2", "u3", "u4"]}
          />
        </div>
      </div>
    </>
  );
};

export default SimulationPage;
