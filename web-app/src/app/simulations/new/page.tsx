"use client"
import { useState } from 'react';
import { useRouter } from 'next/navigation';

export default function NewSimulationPage() {
  const router = useRouter();
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  
  const [formData, setFormData] = useState({
    tle_line1: '',
    tle_line2: '',
    start_date_time: '',
    end_date_time: '',
    control_time_step: '500',
    ground_station_lat: '52.515133',
    ground_station_lon: '13.323456',
    ground_station_alt: '50.0',
    ground_station_elevation: '95.0',
  });

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setLoading(true);
    setError(null);

    try {
      const response = await fetch('/api/v1/simulations', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          tle_line1: formData.tle_line1.trim(),
          tle_line2: formData.tle_line2.trim(),
          start_date_time: formData.start_date_time,
          end_date_time: formData.end_date_time,
          control_time_step: parseInt(formData.control_time_step),
          ground_station_lla: {
            lat: parseFloat(formData.ground_station_lat),
            lon: parseFloat(formData.ground_station_lon),
            alt: parseFloat(formData.ground_station_alt),
          },
          ground_station_elevation: parseFloat(formData.ground_station_elevation),
        }),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.error || 'Failed to create simulation');
      }

      const data = await response.json();
      // Redirect to the simulation list page
      router.push('/');
    } catch (err) {
      setError(err instanceof Error ? err.message : 'An error occurred');
    } finally {
      setLoading(false);
    }
  };

  const handleChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
    setFormData({
      ...formData,
      [e.target.name]: e.target.value,
    });
  };

  // Set default date/time values (current time and 20 minutes later)
  const now = new Date();
  const later = new Date(now.getTime() + 20 * 60 * 1000); // 20 minutes later
  
  const formatDateTime = (date: Date) => {
    return date.toISOString().slice(0, 16); // Format: YYYY-MM-DDTHH:mm
  };

  return (
    <div className="page-container" style={{ maxWidth: '800px', margin: '0 auto' }}>
      <h1>Create New Simulation</h1>
      
      {error && (
        <div style={{
          padding: '1rem',
          marginBottom: '1rem',
          backgroundColor: '#fee2e2',
          color: '#991b1b',
          borderRadius: '0.5rem',
          border: '1px solid #fecaca'
        }}>
          {error}
        </div>
      )}

      <form onSubmit={handleSubmit} style={{ display: 'flex', flexDirection: 'column', gap: '1.5rem' }}>
        <div>
          <label style={{ display: 'block', marginBottom: '0.5rem', fontWeight: '500' }}>
            TLE Line 1 *
          </label>
          <textarea
            name="tle_line1"
            value={formData.tle_line1}
            onChange={handleChange}
            required
            rows={2}
            style={{
              width: '100%',
              padding: '0.5rem',
              border: '1px solid #d1d5db',
              borderRadius: '0.375rem',
              fontFamily: 'monospace',
              fontSize: '0.875rem',
              color: '#000000',
            }}
            placeholder="1 44412U 19038AC  23177.36594369  .00027390  00000+0  92824-3 0  9999"
          />
        </div>

        <div>
          <label style={{ display: 'block', marginBottom: '0.5rem', fontWeight: '500' }}>
            TLE Line 2 *
          </label>
          <textarea
            name="tle_line2"
            value={formData.tle_line2}
            onChange={handleChange}
            required
            rows={2}
            style={{
              width: '100%',
              padding: '0.5rem',
              border: '1px solid #d1d5db',
              borderRadius: '0.375rem',
              fontFamily: 'monospace',
              fontSize: '0.875rem',
              color: '#000000',
            }}
            placeholder="2 44412  97.6739 158.3078 0012642 270.3777  89.6015 15.30455419219567"
          />
        </div>

        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '1rem' }}>
          <div>
            <label style={{ display: 'block', marginBottom: '0.5rem', fontWeight: '500' }}>
              Start Date/Time *
            </label>
            <input
              type="datetime-local"
              name="start_date_time"
              value={formData.start_date_time || formatDateTime(now)}
              onChange={handleChange}
              required
              style={{
                width: '100%',
                padding: '0.5rem',
                border: '1px solid #d1d5db',
                borderRadius: '0.375rem',
                color: '#000000',
              }}
            />
          </div>

          <div>
            <label style={{ display: 'block', marginBottom: '0.5rem', fontWeight: '500' }}>
              End Date/Time *
            </label>
            <input
              type="datetime-local"
              name="end_date_time"
              value={formData.end_date_time || formatDateTime(later)}
              onChange={handleChange}
              required
              style={{
                width: '100%',
                padding: '0.5rem',
                border: '1px solid #d1d5db',
                borderRadius: '0.375rem',
                color: '#000000',
              }}
            />
          </div>
        </div>

        <div>
          <label style={{ display: 'block', marginBottom: '0.5rem', fontWeight: '500' }}>
            Control Time Step (milliseconds) *
          </label>
          <input
            type="number"
            name="control_time_step"
            value={formData.control_time_step}
            onChange={handleChange}
            required
            min="1"
            style={{
              width: '100%',
              padding: '0.5rem',
              border: '1px solid #d1d5db',
              borderRadius: '0.375rem',
              color: '#000000',
            }}
          />
        </div>

        <div style={{ borderTop: '1px solid #e5e7eb', paddingTop: '1.5rem' }}>
          <h2 style={{ fontSize: '1.25rem', marginBottom: '1rem' }}>Ground Station</h2>
          
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr 1fr', gap: '1rem', marginBottom: '1rem' }}>
            <div>
              <label style={{ display: 'block', marginBottom: '0.5rem', fontWeight: '500' }}>
                Latitude (degrees) *
              </label>
              <input
                type="number"
                name="ground_station_lat"
                value={formData.ground_station_lat}
                onChange={handleChange}
                required
                step="any"
                style={{
                  width: '100%',
                  padding: '0.5rem',
                  border: '1px solid #d1d5db',
                  borderRadius: '0.375rem',
                  color: '#000000',
                }}
              />
            </div>

            <div>
              <label style={{ display: 'block', marginBottom: '0.5rem', fontWeight: '500' }}>
                Longitude (degrees) *
              </label>
              <input
                type="number"
                name="ground_station_lon"
                value={formData.ground_station_lon}
                onChange={handleChange}
                required
                step="any"
                style={{
                  width: '100%',
                  padding: '0.5rem',
                  border: '1px solid #d1d5db',
                  borderRadius: '0.375rem',
                  color: '#000000',
                }}
              />
            </div>

            <div>
              <label style={{ display: 'block', marginBottom: '0.5rem', fontWeight: '500' }}>
                Altitude (meters) *
              </label>
              <input
                type="number"
                name="ground_station_alt"
                value={formData.ground_station_alt}
                onChange={handleChange}
                required
                step="any"
                style={{
                  width: '100%',
                  padding: '0.5rem',
                  border: '1px solid #d1d5db',
                  borderRadius: '0.375rem',
                  color: '#000000',
                }}
              />
            </div>
          </div>

          <div>
            <label style={{ display: 'block', marginBottom: '0.5rem', fontWeight: '500' }}>
              Elevation Angle (degrees) *
            </label>
            <input
              type="number"
              name="ground_station_elevation"
              value={formData.ground_station_elevation}
              onChange={handleChange}
              required
              step="any"
              min="0"
              max="180"
              style={{
                width: '100%',
                padding: '0.5rem',
                border: '1px solid #d1d5db',
                borderRadius: '0.375rem',
                color: '#000000',
              }}
            />
          </div>
        </div>

        <div style={{ display: 'flex', gap: '1rem', justifyContent: 'flex-end', marginTop: '1rem' }}>
          <button
            type="button"
            onClick={() => router.back()}
            style={{
              padding: '0.5rem 1.5rem',
              border: '1px solid #d1d5db',
              borderRadius: '0.375rem',
              backgroundColor: 'white',
              cursor: 'pointer',
              color: '#000000',
            }}
          >
            Cancel
          </button>
          <button
            type="submit"
            disabled={loading}
            style={{
              padding: '0.5rem 1.5rem',
              backgroundColor: loading ? '#9ca3af' : '#0070f3',
              color: 'white',
              border: 'none',
              borderRadius: '0.375rem',
              cursor: loading ? 'not-allowed' : 'pointer',
              fontWeight: '500'
            }}
          >
            {loading ? 'Creating...' : 'Create Simulation'}
          </button>
        </div>
      </form>
    </div>
  );
}

