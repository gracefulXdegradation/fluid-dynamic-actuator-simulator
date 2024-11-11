import { readdirSync } from 'fs'
import path from 'path';

export async function GET() {
  try {
    const dataDir = path.join('/data/fds');
    const entries = readdirSync(dataDir, { withFileTypes: true });

    // Filter for directories only and map them to their names
    const subfolderNames = entries
      .filter((entry) => entry.isDirectory())
      .map((folder) => folder.name);

    return Response.json({ entries: subfolderNames });
  } catch (error) {
    // Handle errors, such as if the directory doesn't exist
    console.error(error);
    return new Response(`Failed to retrieve the entries`, {
      status: 500,
    })
  }
}