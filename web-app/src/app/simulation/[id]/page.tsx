"use client"
import { useState, useEffect } from 'react';
import { useParams } from 'next/navigation'
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';

interface Data {
  i_r: number[][]; // Array of 3 arrays for x, y, z components of radius vector
  i_v: number[][]; // Array of 3 arrays for x, y, z components of velocity vector
  t: number[]; // Array of timestamps
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
  const { r, v, t } = data;

  const radiusData = [
    { name: 'X', data: r[0] },
    { name: 'Y', data: r[1] },
    { name: 'Z', data: r[2] },
  ];

  const velocityData = [
    { name: 'X', data: v[0] },
    { name: 'Y', data: v[1] },
    { name: 'Z', data: v[2] },
  ];

  const formatDataForChart = (vectorData: { name: string; data: number[] }[]) => {
    return vectorData[0].data.map((_, index) => {
      // const obj: any = { timestamp: t[index] };
      const obj = { timestamp: new Date().getTime() };
      vectorData.forEach((vector) => {
        obj[vector.name] = vector.data[index];
      });
      return obj;
    });
  };

  const radiusChartData = formatDataForChart(radiusData);
  const velocityChartData = formatDataForChart(velocityData);

  return (
    <div>
      <h1>Simulation {id}</h1>
      <div>
        <h2>Radius Vector</h2>
        <ResponsiveContainer width="100%" height={400}>
          <LineChart data={radiusChartData}>
            <CartesianGrid strokeDasharray="3 3" />
            <XAxis dataKey="timestamp" />
            <YAxis />
            <Tooltip />
            <Legend />
            {radiusData.map((data) => (
              <Line key={data.name} type="monotone" dataKey={data.name} stroke="#8884d8" />
            ))}
          </LineChart>
        </ResponsiveContainer>
      </div>
      <div>
        <h2>Velocity Vector</h2>
        <ResponsiveContainer width="100%" height={400}>
          <LineChart data={velocityChartData}>
            <CartesianGrid strokeDasharray="3 3" />
            <XAxis dataKey="timestamp" />
            <YAxis />
            <Tooltip />
            <Legend />
            {velocityData.map((data) => (
              <Line key={data.name} type="monotone" dataKey={data.name} stroke="#82ca9d" />
            ))}
          </LineChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
};

export default SimulationPage;
