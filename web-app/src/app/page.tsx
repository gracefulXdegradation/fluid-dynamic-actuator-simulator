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
    <div className="page-container">
      <h1>Simulations</h1>
      {loading ? (
        <div className="loading-indicator"/>
      ) : entries.length > 0 ? (
        <ul>
          {entries.map((entry) => (
            <li key={entry}>
              <Link href={`/simulation/${entry}`}>
                #{entry}
              </Link>
            </li>
          ))}
        </ul>
      ) : (
        <div>No entries available.</div>
      )}
    </div>
  );
};

export default SimulationsPage;
