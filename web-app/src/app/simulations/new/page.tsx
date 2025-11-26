"use client"
import { useState } from 'react';
import { useRouter } from 'next/navigation';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Textarea } from '@/components/ui/textarea';
import { Label } from '@/components/ui/label';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';

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
    <div className="page-container">
      <div className="w-full max-w-3xl mx-auto">
        <h1 className="text-2xl font-bold mb-6">Create New Simulation</h1>
        
        {error && (
          <Alert variant="destructive" className="mb-6">
            <AlertDescription>{error}</AlertDescription>
          </Alert>
        )}

        <Card>
          <CardHeader>
            <CardTitle>Simulation Details</CardTitle>
            <CardDescription>Enter the parameters for your new simulation</CardDescription>
          </CardHeader>
          <CardContent>
            <form onSubmit={handleSubmit} className="space-y-6">
              <div className="space-y-2">
                <Label htmlFor="tle_line1">
                  TLE Line 1 *
                </Label>
                <Textarea
                  id="tle_line1"
                  name="tle_line1"
                  value={formData.tle_line1}
                  onChange={handleChange}
                  required
                  rows={1}
                  className="font-mono text-sm"
                  placeholder="1 44412U 19038AC  23177.36594369  .00027390  00000+0  92824-3 0  9999"
                />
              </div>

              <div className="space-y-2">
                <Label htmlFor="tle_line2">
                  TLE Line 2 *
                </Label>
                <Textarea
                  id="tle_line2"
                  name="tle_line2"
                  value={formData.tle_line2}
                  onChange={handleChange}
                  required
                  rows={1}
                  className="font-mono text-sm"
                  placeholder="2 44412  97.6739 158.3078 0012642 270.3777  89.6015 15.30455419219567"
                />
              </div>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <div className="space-y-2">
                  <Label htmlFor="start_date_time">
                    Start Date/Time *
                  </Label>
                  <Input
                    id="start_date_time"
                    type="datetime-local"
                    name="start_date_time"
                    value={formData.start_date_time || formatDateTime(now)}
                    onChange={handleChange}
                    required
                  />
                </div>

                <div className="space-y-2">
                  <Label htmlFor="end_date_time">
                    End Date/Time *
                  </Label>
                  <Input
                    id="end_date_time"
                    type="datetime-local"
                    name="end_date_time"
                    value={formData.end_date_time || formatDateTime(later)}
                    onChange={handleChange}
                    required
                  />
                </div>
              </div>

              <div className="space-y-2">
                <Label htmlFor="control_time_step">
                  Control Time Step (milliseconds) *
                </Label>
                <Input
                  id="control_time_step"
                  type="number"
                  name="control_time_step"
                  value={formData.control_time_step}
                  onChange={handleChange}
                  required
                  min="1"
                />
              </div>

              <div className="border-t pt-6 space-y-6">
                <div>
                  <h2 className="text-xl font-semibold mb-4">Ground Station</h2>
                  
                  <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
                    <div className="space-y-2">
                      <Label htmlFor="ground_station_lat">
                        Latitude (degrees) *
                      </Label>
                      <Input
                        id="ground_station_lat"
                        type="number"
                        name="ground_station_lat"
                        value={formData.ground_station_lat}
                        onChange={handleChange}
                        required
                        step="any"
                      />
                    </div>

                    <div className="space-y-2">
                      <Label htmlFor="ground_station_lon">
                        Longitude (degrees) *
                      </Label>
                      <Input
                        id="ground_station_lon"
                        type="number"
                        name="ground_station_lon"
                        value={formData.ground_station_lon}
                        onChange={handleChange}
                        required
                        step="any"
                      />
                    </div>

                    <div className="space-y-2">
                      <Label htmlFor="ground_station_alt">
                        Altitude (meters) *
                      </Label>
                      <Input
                        id="ground_station_alt"
                        type="number"
                        name="ground_station_alt"
                        value={formData.ground_station_alt}
                        onChange={handleChange}
                        required
                        step="any"
                      />
                    </div>
                  </div>

                  <div className="space-y-2">
                    <Label htmlFor="ground_station_elevation">
                      Elevation Angle (degrees) *
                    </Label>
                    <Input
                      id="ground_station_elevation"
                      type="number"
                      name="ground_station_elevation"
                      value={formData.ground_station_elevation}
                      onChange={handleChange}
                      required
                      step="any"
                      min="0"
                      max="180"
                    />
                  </div>
                </div>
              </div>

              <div className="flex gap-4 justify-end pt-4">
                <Button
                  type="button"
                  variant="outline"
                  onClick={() => router.back()}
                >
                  Cancel
                </Button>
                <Button
                  type="submit"
                  disabled={loading}
                >
                  {loading ? 'Creating...' : 'Create Simulation'}
                </Button>
              </div>
            </form>
          </CardContent>
        </Card>
      </div>
    </div>
  );
}

