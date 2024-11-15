"use client"
import { useState, useEffect } from 'react';
import { useParams } from 'next/navigation'
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';

const colors = ["#8884d8", "#82ca9d", "#ffa600", "#ff6361"];

interface Data {
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
  const { t, a_control_torque, d } = data;

  const ts = t[0];

  const formatDataForChart = (vectorData: { name: string; data: number[] }[], timestamps: number[]) => {
    const data = timestamps.reduce((acc, ti, index) => {
      const obj = {
        name: new Date(ti * 1000),
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

  return (
    <div>
      <h1>Simulation {id}</h1>
      <div>
        <h2>Required control torque</h2>
        <ResponsiveContainer width="50%" height={400}>
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
    </div>
  );
};

export default SimulationPage;
