import { promises as fs } from 'fs';
import path from 'path';
import { parse } from 'csv-parse';

const parseCSV = async (filePath: string) => {
      // Read the CSV file
      const fileContent = await fs.readFile(filePath, 'utf-8');

      // Parse CSV data into JSON
      let rows: number[][];
      parse(fileContent, { cast: true }, (err, records: number[][]) => {
        if (err) {
          throw err;
        }
        
        if (records[0] && records[0].length === 1) {
          rows = records.map(r => r[0]);
        } else {
          rows = records;
        }
      });
  
      // Wait for parsing to complete
      await new Promise(resolve => setTimeout(resolve, 100)); // minor delay to complete parsing

      return rows;
}

export async function GET(request: Request, { params }: { params: Promise<{ timestamp: string; }> }) {
  const {timestamp} = await params;

  try {
    const euler_angles = await parseCSV(path.join('/data/fds', timestamp, `euler_angles.csv`));
    const a_control_torque = await parseCSV(path.join('/data/fds', timestamp, `a_control_torque.csv`));
    const ang_mom_body_frame = await parseCSV(path.join('/data/fds', timestamp, `ang_mom_body_frame.csv`));
    const a_command = await parseCSV(path.join('/data/fds', timestamp, `a_command.csv`));
    const state = await parseCSV(path.join('/data/fds', timestamp, `state.csv`));
    const t = await parseCSV(path.join('/data/fds', timestamp, `t.csv`));
    const distance = await parseCSV(path.join('/data/fds', timestamp, `distance.csv`));

    return Response.json({
      euler_angles,
      ang_mom_body_frame,
      a_control_torque,
      a_command,
      state,
      t,
      d: distance});
  } catch (error) {
    console.error('Error reading file:', error);
    return new Response(`Path does not exist: ${timestamp}/${filename}`, {
      status: 500,
    })
  }
}