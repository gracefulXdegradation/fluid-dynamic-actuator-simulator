"use client"
import { useState, useEffect } from 'react';
import { useParams } from 'next/navigation'
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';

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
  const { t, d } = data;

  const ts = t[0];

  // const radiusData = [
  //   { name: 'X', data: r[0] },
  //   { name: 'Y', data: r[1] },
  //   { name: 'Z', data: r[2] },
  // ];

  // const velocityData = [
  //   { name: 'X', data: v[0] },
  //   { name: 'Y', data: v[1] },
  //   { name: 'Z', data: v[2] },
  // ];

  const formatDataForChart = (vectorData: { name: string; data: number[] }[], timestamps: number[]) => {
    return timestamps.reduce((acc, ti, index) => {
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

    // return vectorData[0].data.map((_, index) => {
    //   const obj = { timestamp: new Date(timestamps[index] * 1000) };
    //   vectorData.forEach((vector) => {
    //     obj[vector.name] = vector.data[index];
    //   });
    //   return obj;
    // });
  };

  // const radiusChartData = formatDataForChart(radiusData);
  // const velocityChartData = formatDataForChart(velocityData);
  const distanceChartData = formatDataForChart([{
    name: 'Distance', data: d[0]
  }], ts);

  return (
    <div>
      <h1>Simulation {id}</h1>
      {/* <div>
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
      </div> */}
      <div>
        <h2>Distance</h2>
        <ResponsiveContainer width="100%" height={400}>
          <LineChart data={distanceChartData}>
            <CartesianGrid strokeDasharray="3 3" />
            <XAxis dataKey="timestamp" />
            <YAxis />
            <Tooltip />
            <Legend />
            {/* {distanceChartData.map((data) => ( */}
              <Line type="monotone" dataKey="Distance" stroke="#82ca9d" />
            {/* ))} */}
          </LineChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
};

export default SimulationPage;
