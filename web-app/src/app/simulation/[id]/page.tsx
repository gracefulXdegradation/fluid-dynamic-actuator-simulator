"use client"
import { useState, useEffect } from 'react';
import { useParams } from 'next/navigation'
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';
import LineGraph from '@/components/LineGraph';

const colors = ["#8884d8", "#82ca9d", "#ffa600", "#ff6361"];

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
  const angularMomentum = [7, 9, 11 ,13].map(i => state[i]);
  const ts = t[0];
  
  // window.angularRate = angularRate;
  // window.angularRateAbs = angularRateAbs;

  const formatDataForChart = (vectorData: { name: string; data: number[] }[], timestamps: number[]) => {
    const data = timestamps.reduce((acc, ti, index) => {
      const obj = {
        name: new Date(ti),
        t: ti        
      };

      vectorData.forEach((vector) => {
        obj[vector.name] = vector.data[index];
      });

      acc.push(obj)

      return acc;
    }, [] as {name: Date; t: number; [key: string]: number}[])

    return {
      data,
      meta: vectorData.map((v, i) => ({
        key: v.name,
        stroke: v.stroke || colors[i]
      }))
    };
  };

  const distanceChartData = formatDataForChart([{
    name: 'Distance',
    data: d[0]
  }], ts);

  const controlTorqueData = formatDataForChart([{
    name: 'Tau1',
    data: a_control_torque[0]
  }, {
    name: 'Tau2',
    data: a_control_torque[1]
  }, {
    name: 'Tau3',
    data: a_control_torque[2]
  }, {
    name: 'Tau4',
    data: a_control_torque[3]
  }], ts);

  const commandTorqueData = formatDataForChart([{
    name: 'u1',
    data: a_command[0]
  }, {
    name: 'u2',
    data: a_command[1]
  }, {
    name: 'u3',
    data: a_command[2]
  }, {
    name: 'u4',
    data: a_command[3]
  }], ts);

  const angularRateData = formatDataForChart([{
    name: 'Omega x',
    data: angularRate[0]
  }, {
    name: 'Omega y',
    data: angularRate[1]
  }, {
    name: 'Omega z',
    data: angularRate[2]
  }, {
    name: '|Omega|',
    data: angularRateAbs
  }], ts);

  const angularMomentumData = formatDataForChart([{
    name: 'h1',
    data: angularMomentum[0]
  }, {
    name: 'h2',
    data: angularMomentum[1]
  }, {
    name: 'h3',
    data: angularMomentum[2]
  }, {
    name: 'h4',
    data: angularMomentum[3]
  }], ts);

  console.log('angularMomentumData', angularMomentumData)

  const angularMomentumBodyFrameData = formatDataForChart([{
    name: 'h_x',
    data: ang_mom_body_frame[0]
  }, {
    name: 'h_y',
    data: ang_mom_body_frame[1]
  }, {
    name: 'h_z',
    data: ang_mom_body_frame[2]
  }], ts);

  const eulerAnglesData = formatDataForChart([{
    name: 'y',
    // clamp the value between 0 and 0.2 degrees
    data: euler_angles[1].map(rad => Math.min(Math.max(rad2deg(rad), 0), 0.2))
  }], ts);

  return (
    <div>
      <h1>Simulation {id}</h1>
      <div>
      <LineGraph
        timestamps={ts}
        values={ang_mom_body_frame}
        graphNames={["h_x", "h_y", "h_z"]}
      />
        <h2>Required control torque</h2>
        {/* <ResponsiveContainer width="50%" height={400}>
          <LineChart data={controlTorqueData.data}>
            <CartesianGrid strokeDasharray="3 3" opacity={0.5} />
            <XAxis dataKey="timestamp" />
            <YAxis />
            <Tooltip />
            <Legend />
            {controlTorqueData.meta.map(({key, stroke}) => (
              <Line type="monotone" key={key} dataKey={key} stroke={stroke} />
            ))}
          </LineChart>
        </ResponsiveContainer>
      </div>
      <div>
        <h2>Body angular rate w.r.t. body frame</h2>
        <ResponsiveContainer width="50%" height={400}>
          <LineChart data={angularRateData.data}>
            <CartesianGrid strokeDasharray="3 3" opacity={0.5} />
            <XAxis dataKey="timestamp" />
            <YAxis />
            <Tooltip />
            <Legend />
            {angularRateData.meta.map(({key, stroke}) => (
              <Line type="monotone" key={key} dataKey={key} stroke={stroke} />
            ))}
          </LineChart>
        </ResponsiveContainer>
      </div>
      <div>
        <h2>Actuator angular momentum in actuator frame</h2>
        <ResponsiveContainer width="50%" height={400}>
          <LineChart data={angularMomentumData.data}>
            <CartesianGrid strokeDasharray="3 3" opacity={0.5} />
            <XAxis dataKey="timestamp" />
            <YAxis />
            <Tooltip />
            <Legend />
            {angularMomentumData.meta.map(({key, stroke}) => (
              <Line type="monotone" key={key} dataKey={key} stroke={stroke} />
            ))}
          </LineChart>
        </ResponsiveContainer>
      </div>
      <div>
        <h2>Actuator angular momentum in body frame</h2>
        <ResponsiveContainer width="50%" height={400}>
          <LineChart data={angularMomentumBodyFrameData.data}>
            <CartesianGrid strokeDasharray="3 3" opacity={0.5} />
            <XAxis dataKey="timestamp" />
            <YAxis />
            <Tooltip />
            <Legend />
            {angularMomentumBodyFrameData.meta.map(({key, stroke}) => (
              <Line type="monotone" key={key} dataKey={key} stroke={stroke} />
            ))}
          </LineChart>
        </ResponsiveContainer>
      </div>
      <div>
        <h2>Distance</h2>
        <ResponsiveContainer width="50%" height={400}>
          <LineChart data={distanceChartData.data}>
            <CartesianGrid strokeDasharray="5 5" opacity={0.5} />
            <XAxis dataKey="timestamp" />
            <YAxis />
            <Tooltip />
            <Legend />
            {distanceChartData.meta.map(({key, stroke}) => (
              <Line type="monotone" key={key} dataKey={key} stroke={stroke} />
            ))}
          </LineChart>
        </ResponsiveContainer>
      </div>
      <div>
        <h2>Attitude error angle</h2>
        <ResponsiveContainer width="50%" height={400}>
          <LineChart data={eulerAnglesData.data}>
            <CartesianGrid strokeDasharray="5 5" opacity={0.5} />
            <XAxis dataKey="timestamp" />
            <YAxis />
            <Tooltip />
            <Legend />
            {eulerAnglesData.meta.map(({key, stroke}) => (
              <Line type="monotone" key={key} dataKey={key} stroke={stroke} />
            ))}
          </LineChart>
        </ResponsiveContainer>
      </div>
      <div>
        <h2>Actuator commands</h2>
        <ResponsiveContainer width="50%" height={400}>
          <LineChart data={commandTorqueData.data}>
            <CartesianGrid strokeDasharray="5 5" opacity={0.5} />
            <XAxis dataKey="timestamp" />
            <YAxis />
            <Tooltip />
            <Legend />
            {commandTorqueData.meta.map(({key, stroke}) => (
              <Line type="monotone" key={key} dataKey={key} stroke={stroke} />
            ))}
          </LineChart>
        </ResponsiveContainer> */}
      </div>
    </div>
  );
};

export default SimulationPage;
