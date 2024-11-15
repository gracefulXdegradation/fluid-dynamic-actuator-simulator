"use client"
import { useState, useEffect } from 'react';
import Link from 'next/link';

const SimulationsPage = () => {
  const [entries, setEntries] = useState<string[]>([]);
  const [loading, setLoading] = useState(true);
  

  // Fetch entries from the API
  const fetchEntries = async () => {
    try {
      const response = await fetch('/api/v1/entries');
      const data = await response.json();
      if (Array.isArray(data.entries)) {
        setEntries(data.entries);
      } else {
        setEntries([]);
      }
    } catch (error) {
      console.error('Error fetching entries:', error);
      setEntries([]);
    } finally {
      setLoading(false);
    }
  };

  // Set up polling with a 5-second interval
  useEffect(() => {
    fetchEntries(); // Initial fetch
    const interval = setInterval(fetchEntries, 5000);
    return () => clearInterval(interval); // Clean up on unmount
  }, []);

  return (
    <div style={{ textAlign: 'center' }}>
      <h1>Simulation Entries</h1>
      {loading ? (
        <div className="loading-indicator">Loading...</div>
      ) : entries.length > 0 ? (
        <ul>
          {entries.map((entry) => (
            <li key={entry}>
              <Link href={`/simulation/${entry}`}>
                {entry}
              </Link>
            </li>
          ))}
        </ul>
      ) : (
        <div>No entries available.</div>
      )}
      <style jsx>{`
        .loading-indicator {
          border: 4px solid #f3f3f3;
          border-top: 4px solid #3498db;
          border-radius: 50%;
          width: 30px;
          height: 30px;
          animation: spin 1s linear infinite;
          margin: 20px auto;
        }
        @keyframes spin {
          0% {
            transform: rotate(0deg);
          }
          100% {
            transform: rotate(360deg);
          }
        }
      `}</style>
    </div>
  );
};

export default SimulationsPage;
