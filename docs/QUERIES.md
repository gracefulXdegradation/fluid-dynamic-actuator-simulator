# Database Inspection Commands

## Quick Commands

### Connect to Database

```bash
# Find the postgres container name first
docker ps | grep postgres

# Connect (replace 'fds-postgres-1' with your actual container name)
docker exec -it fds-postgres-1 psql -U fds_user -d fds_db
```

Or if you have `psql` installed locally:

```bash
psql postgresql://fds_user:fds_password@localhost:5432/fds_db
```

### List All Tables

```sql
\dt
```

### View Table Structure

```sql
-- View simulations table structure
\d simulations

-- View simulation_metrics table structure
\d simulation_metrics
```

### List Custom Types (Enums)

```sql
-- List all custom types
\dT

-- View simulation_status enum details
\dT+ simulation_status
```

### View Sample Data

```sql
-- View all simulations
SELECT * FROM simulations ORDER BY created_at DESC LIMIT 10;

-- View simulation with its metrics count
SELECT 
    s.id,
    s.status,
    s.created_at,
    s.started_at,
    s.completed_at,
    COUNT(sm.id) as metric_count
FROM simulations s
LEFT JOIN simulation_metrics sm ON s.id = sm.simulation_id
GROUP BY s.id
ORDER BY s.created_at DESC
LIMIT 10;

-- View metrics for a specific simulation
SELECT * FROM simulation_metrics 
WHERE simulation_id = 'YOUR_SIMULATION_ID'
ORDER BY step_index
LIMIT 10;
```

### Get Table Sizes

```sql
-- Table sizes
SELECT 
    schemaname,
    tablename,
    pg_size_pretty(pg_total_relation_size(schemaname||'.'||tablename)) AS size
FROM pg_tables
WHERE schemaname = 'public'
ORDER BY pg_total_relation_size(schemaname||'.'||tablename) DESC;
```

### View Indexes

```sql
-- List all indexes
\di

-- Or query
SELECT 
    tablename,
    indexname,
    indexdef
FROM pg_indexes
WHERE schemaname = 'public'
ORDER BY tablename, indexname;
```

### Export Schema

```sql
-- Export just the schema (no data)
\dn+  -- List schemas
\dt+  -- List tables with sizes
```

## One-Line Commands (without entering psql)

```bash
# List all tables
docker exec fds-postgres-1 psql -U fds_user -d fds_db -c "\dt"

# Describe simulations table
docker exec fds-postgres-1 psql -U fds_user -d fds_db -c "\d simulations"

# Describe simulation_metrics table
docker exec fds-postgres-1 psql -U fds_user -d fds_db -c "\d simulation_metrics"

# Count simulations by status
docker exec fds-postgres-1 psql -U fds_user -d fds_db -c "SELECT status, COUNT(*) FROM simulations GROUP BY status;"

# View recent simulations
docker exec fds-postgres-1 psql -U fds_user -d fds_db -c "SELECT id, status, created_at FROM simulations ORDER BY created_at DESC LIMIT 5;"

# View the latest time stamps
docker exec fds-postgres-1 psql -U fds_user -d fds_db -c "SELECT id, simulation_id, step_index, timestamp, euler_angles FROM simulation_metrics ORDER BY step_index DESC LIMIT 10;"
```

