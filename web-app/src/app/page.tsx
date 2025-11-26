"use client"
import { useState, useEffect } from 'react';
import Link from 'next/link';

type SimulationStatus = 'scheduled' | 'running' | 'completed' | 'failed';

interface Simulation {
  id: string;
  status: SimulationStatus;
  created_at: string;
  started_at: string | null;
  completed_at: string | null;
  error_message: string | null;
}

const SimulationsPage = () => {
  const [simulations, setSimulations] = useState<Simulation[]>([]);
  const [loading, setLoading] = useState(true);

  // Fetch simulations from the API
  const fetchSimulations = async () => {
    try {
      const response = await fetch('/api/v1/simulations');
      const data = await response.json();
      if (Array.isArray(data.simulations)) {
        setSimulations(data.simulations);
      } else {
        setSimulations([]);
      }
    } catch (error) {
      console.error('Error fetching simulations:', error);
      setSimulations([]);
    } finally {
      setLoading(false);
    }
  };

  // Set up polling with a 5-second interval
  useEffect(() => {
    fetchSimulations(); // Initial fetch
    const interval = setInterval(fetchSimulations, 5000);
    return () => clearInterval(interval); // Clean up on unmount
  }, []);

  const getStatusColor = (status: SimulationStatus) => {
    switch (status) {
      case 'scheduled':
        return 'text-blue-600 bg-blue-100';
      case 'running':
        return 'text-yellow-600 bg-yellow-100';
      case 'completed':
        return 'text-green-600 bg-green-100';
      case 'failed':
        return 'text-red-600 bg-red-100';
      default:
        return 'text-gray-600 bg-gray-100';
    }
  };

  const formatDate = (dateString: string | null) => {
    if (!dateString) return 'N/A';
    return new Date(dateString).toLocaleString();
  };

  return (
    <div className="page-container">
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '2rem' }}>
        <h1>Simulations</h1>
        <Link 
          href="/simulations/new"
          style={{
            padding: '0.5rem 1rem',
            backgroundColor: '#0070f3',
            color: 'white',
            borderRadius: '0.5rem',
            textDecoration: 'none',
            fontWeight: '500'
          }}
        >
          + New Simulation
        </Link>
      </div>
      
      {loading ? (
        <div className="loading-indicator"/>
      ) : simulations.length > 0 ? (
        <div style={{ overflowX: 'auto' }}>
          <table style={{ width: '100%', borderCollapse: 'collapse' }}>
            <thead>
              <tr style={{ borderBottom: '2px solid #e5e7eb' }}>
                <th style={{ padding: '0.75rem', textAlign: 'left' }}>ID</th>
                <th style={{ padding: '0.75rem', textAlign: 'left' }}>Status</th>
                <th style={{ padding: '0.75rem', textAlign: 'left' }}>Created</th>
                <th style={{ padding: '0.75rem', textAlign: 'left' }}>Started</th>
                <th style={{ padding: '0.75rem', textAlign: 'left' }}>Completed</th>
                <th style={{ padding: '0.75rem', textAlign: 'left' }}>Actions</th>
              </tr>
            </thead>
            <tbody>
              {simulations.map((sim) => (
                <tr key={sim.id} style={{ borderBottom: '1px solid #e5e7eb' }}>
                  <td style={{ padding: '0.75rem', fontFamily: 'monospace', fontSize: '0.875rem' }}>
                    {sim.id.substring(0, 8)}...
                  </td>
                  <td style={{ padding: '0.75rem' }}>
                    <span 
                      style={{
                        padding: '0.25rem 0.75rem',
                        borderRadius: '9999px',
                        fontSize: '0.875rem',
                        fontWeight: '500',
                        display: 'inline-block'
                      }}
                      className={getStatusColor(sim.status)}
                    >
                      {sim.status}
                    </span>
                  </td>
                  <td style={{ padding: '0.75rem', fontSize: '0.875rem' }}>
                    {formatDate(sim.created_at)}
                  </td>
                  <td style={{ padding: '0.75rem', fontSize: '0.875rem' }}>
                    {formatDate(sim.started_at)}
                  </td>
                  <td style={{ padding: '0.75rem', fontSize: '0.875rem' }}>
                    {formatDate(sim.completed_at)}
                  </td>
                  <td style={{ padding: '0.75rem' }}>
                    <Link 
                      href={`/simulations/${sim.id}`}
                      style={{
                        color: '#0070f3',
                        textDecoration: 'none',
                        fontWeight: '500'
                      }}
                    >
                      View
                    </Link>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      ) : (
        <div style={{ textAlign: 'center', padding: '3rem', color: '#6b7280' }}>
          <p>No simulations available.</p>
          <Link 
            href="/simulations/new"
            style={{
              display: 'inline-block',
              marginTop: '1rem',
              padding: '0.5rem 1rem',
              backgroundColor: '#0070f3',
              color: 'white',
              borderRadius: '0.5rem',
              textDecoration: 'none'
            }}
          >
            Create your first simulation
          </Link>
        </div>
      )}
    </div>
  );
};

export default SimulationsPage;
