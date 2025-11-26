"use client"
import { useState, useEffect } from 'react';
import Link from 'next/link';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Skeleton } from '@/components/ui/skeleton';

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

  const getStatusBadgeProps = (status: SimulationStatus) => {
    switch (status) {
      case 'scheduled':
        return { variant: 'secondary' as const, className: 'bg-blue-100 text-blue-600 hover:bg-blue-100' };
      case 'running':
        return { variant: 'default' as const, className: 'bg-yellow-100 text-yellow-600 hover:bg-yellow-100' };
      case 'completed':
        return { variant: 'default' as const, className: 'bg-green-100 text-green-600 hover:bg-green-100' };
      case 'failed':
        return { variant: 'destructive' as const };
      default:
        return { variant: 'outline' as const };
    }
  };

  const formatDate = (dateString: string | null) => {
    if (!dateString) return 'N/A';
    return new Date(dateString).toLocaleString();
  };

  return (
    <div className="page-container">
      <div className="flex justify-between items-center mb-8 w-full max-w-7xl">
        <h1 className="text-2xl font-bold">Simulations</h1>
        <Button asChild>
          <Link href="/simulations/new">
            + New Simulation
          </Link>
        </Button>
      </div>
      
      {loading ? (
        <Card className="w-full max-w-7xl">
          <CardHeader>
            <Skeleton className="h-8 w-48" />
          </CardHeader>
          <CardContent>
            <div className="space-y-4">
              {[1, 2, 3, 4, 5].map((i) => (
                <Skeleton key={i} className="h-16 w-full" />
              ))}
            </div>
          </CardContent>
        </Card>
      ) : simulations.length > 0 ? (
        <Card className="w-full max-w-7xl">
          <CardHeader>
            <CardTitle>Simulations</CardTitle>
            <CardDescription>View and manage your simulations</CardDescription>
          </CardHeader>
          <CardContent>
            <Table>
              <TableHeader>
                <TableRow>
                  <TableHead>ID</TableHead>
                  <TableHead>Status</TableHead>
                  <TableHead>Created</TableHead>
                  <TableHead>Started</TableHead>
                  <TableHead>Completed</TableHead>
                  <TableHead>Actions</TableHead>
                </TableRow>
              </TableHeader>
              <TableBody>
                {simulations.map((sim) => (
                  <TableRow key={sim.id}>
                    <TableCell className="font-mono text-sm">
                      {sim.id.substring(0, 8)}...
                    </TableCell>
                    <TableCell>
                      <Badge {...getStatusBadgeProps(sim.status)}>
                        {sim.status}
                      </Badge>
                    </TableCell>
                    <TableCell className="text-sm">
                      {formatDate(sim.created_at)}
                    </TableCell>
                    <TableCell className="text-sm">
                      {formatDate(sim.started_at)}
                    </TableCell>
                    <TableCell className="text-sm">
                      {formatDate(sim.completed_at)}
                    </TableCell>
                    <TableCell>
                      <Button variant="link" asChild>
                        <Link href={`/simulations/${sim.id}`}>
                          View
                        </Link>
                      </Button>
                    </TableCell>
                  </TableRow>
                ))}
              </TableBody>
            </Table>
          </CardContent>
        </Card>
      ) : (
        <Card className="w-full max-w-7xl">
          <CardContent className="flex flex-col items-center justify-center py-12">
            <p className="text-muted-foreground mb-4">No simulations available.</p>
            <Button asChild>
              <Link href="/simulations/new">
                Create your first simulation
              </Link>
            </Button>
          </CardContent>
        </Card>
      )}
    </div>
  );
};

export default SimulationsPage;
