import { promises as fs } from 'fs';
import path from 'path';
import { parse } from 'csv-parse';

const parseCSV = async (filePath: string) => {
      // Read the CSV file
      const fileContent = await fs.readFile(filePath, 'utf-8');

      // Parse CSV data into JSON
      const rows: number[][] = [];
      parse(fileContent, { cast: true }, (err, records: number[][]) => {
        if (err) {
          throw err;
        }
        rows.push(...records);
      });
  
      // Wait for parsing to complete
      await new Promise(resolve => setTimeout(resolve, 100)); // minor delay to complete parsing

      return rows;
}

export async function GET(request: Request, { params }: { params: Promise<{ timestamp: string; }> }) {
  // const searchParams = request.nextUrl.searchParams
  // const ts = searchParams.get('timestamp')
  // const filename = searchParams.get('filename')
  const {timestamp} = await params;

  try {
    const radiusVectors = await parseCSV(path.join('/data/fds', timestamp, `i_r.csv`));
    const velocityVectors = await parseCSV(path.join('/data/fds', timestamp, `i_v.csv`));

    return Response.json({r: radiusVectors, v: velocityVectors});
  } catch (error) {
    console.error('Error reading file:', error);
    return new Response(`Path does not exist: ${timestamp}/${filename}`, {
      status: 500,
    })
  }
}