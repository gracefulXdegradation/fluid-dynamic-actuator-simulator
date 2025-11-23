// Redis client for the web app
// Uses ioredis for better TypeScript support and async/await

let redisClient: any = null;

export class RedisClient {
  private client: any;

  constructor() {
    if (!redisClient) {
      // Lazy load ioredis only when needed
      const Redis = require('ioredis');
      const redisUrl = process.env.REDIS_URL || 'redis://redis:6379';
      redisClient = new Redis(redisUrl);
    }
    this.client = redisClient;
  }

  async pushJob(simulationId: string): Promise<void> {
    await this.client.lpush('simulation_jobs', simulationId);
  }

  async getQueueLength(): Promise<number> {
    return await this.client.llen('simulation_jobs');
  }
}

