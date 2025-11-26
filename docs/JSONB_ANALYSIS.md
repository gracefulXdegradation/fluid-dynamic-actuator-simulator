# JSONB vs Normalized Columns: Trade-off Analysis for Simulation Parameters

## Executive Summary

**Current Implementation**: Simulation parameters stored as JSONB in `input_parameters` column

**Recommendation**: **JSONB is appropriate for this use case**, but consider adding generated columns for frequently queried fields if search/filter requirements emerge.

---

## Current Data Structure

### JSONB Content (`input_parameters`)
```json
{
  "tle_line1": "1 25544U 98067A   08264.51782528 ...",
  "tle_line2": "2 25544  51.6416 247.4627 0006703 ...",
  "start_date_time": "2023-07-04T14:25:00Z",
  "end_date_time": "2023-07-05T14:25:00Z",
  "control_time_step": 1000,
  "ground_station_lla": {
    "lat": 40.7128,
    "lon": -74.0060,
    "alt": 0.0
  },
  "ground_station_elevation": 5.0
}
```

### Current Query Patterns
- ✅ **Primary key lookup**: `SELECT input_parameters FROM simulations WHERE id = $1`
- ✅ **Full document retrieval**: Always fetch entire JSONB, parse in application
- ❌ **No filtering by parameters**: No queries filter by TLE, dates, or ground station
- ❌ **No aggregations**: No queries count/group by parameter values

---

## Pros of JSONB Approach

### 1. **Schema Flexibility** ✅
- **Easy to add new parameters** without ALTER TABLE migrations
- **Backward compatible**: Old code can ignore new fields
- **Optional fields**: Missing fields can be handled gracefully
- **Nested structures**: Natural representation of `ground_station_lla` object

**Example Benefit**:
```sql
-- Adding new parameter requires no schema change
-- Just update application code to handle new field
{
  "tle_line1": "...",
  "new_parameter": "value"  -- Added without migration
}
```

### 2. **Atomic Updates** ✅
- **Single column update**: Update entire parameter set atomically
- **No partial updates**: Avoids inconsistent state across multiple columns
- **Simpler transactions**: One UPDATE statement instead of multiple

**Current Code**:
```cpp
// Single atomic update
"ON CONFLICT (id) DO UPDATE SET input_parameters = $2::jsonb"
```

### 3. **Application Code Simplicity** ✅
- **Single source of truth**: All parameters in one place
- **Natural JSON serialization**: Easy to serialize/deserialize from API
- **Type validation in application**: Can validate entire structure at once
- **Already using JSON library**: `nlohmann::json` is already in use

### 4. **Storage Efficiency** ✅
- **PostgreSQL compression**: JSONB is stored in binary format, compressed
- **No NULL overhead**: Missing optional fields don't take space
- **Efficient for sparse data**: If parameters vary significantly between simulations

### 5. **Query Capabilities (When Needed)** ✅
- **JSONB operators**: Can query inside JSONB if needed later
  ```sql
  -- Example: Find simulations with specific ground station
  SELECT * FROM simulations 
  WHERE input_parameters->'ground_station_lla'->>'lat' = '40.7128';
  ```
- **GIN indexes**: Can index JSONB fields for fast searches
- **JSONB functions**: Rich set of PostgreSQL JSONB functions available

### 6. **API Integration** ✅
- **Direct API mapping**: JSONB can be directly returned from REST APIs
- **No transformation layer**: No need to convert between DB columns and JSON
- **Versioning**: Can include schema version in JSONB for migration handling

---

## Cons of JSONB Approach

### 1. **Limited Query Performance** ⚠️
- **No native indexes on nested fields**: Can't create B-tree index on `ground_station_lla.lat` directly
- **Full table scans for filtering**: If you need to find simulations by date range, must scan all rows
- **JSONB operator overhead**: Queries like `->` and `->>` are slower than column comparisons

**Impact**: Currently **LOW** - you only query by `id` (primary key)

**Future Risk**: If you need queries like:
```sql
-- This would be slow without proper indexing
SELECT * FROM simulations 
WHERE (input_parameters->>'start_date_time')::timestamp > '2023-01-01';
```

### 2. **No Database-Level Constraints** ⚠️
- **No NOT NULL constraints**: Can't enforce required fields at DB level
- **No CHECK constraints**: Can't validate ranges (e.g., `control_time_step > 0`)
- **No foreign keys**: Can't reference other tables from JSONB fields
- **No data types**: All values stored as JSON types, not native PostgreSQL types

**Current Mitigation**: Application-level validation in `parseInputParameters()`

**Example Risk**:
```json
// This would be accepted by database, but fail at runtime
{
  "control_time_step": -1000,  // Invalid negative value
  "start_date_time": "invalid-date"  // Invalid format
}
```

### 3. **Type Safety** ⚠️
- **Runtime type checking**: Types only validated when parsing JSON
- **No compile-time guarantees**: Can't use PostgreSQL's type system
- **String conversions**: Dates stored as strings, converted in application
- **Precision loss risk**: JSON numbers are double-precision, may lose precision for very large integers

**Current Code**:
```cpp
// Type conversion happens at runtime
params.start_date_time = DateTime::parseDateTime(
    params_json["start_date_time"].get<std::string>()
);
```

### 4. **Migration Complexity** ⚠️
- **Schema evolution**: Need to handle multiple JSONB schema versions in application
- **Data migration**: Changing structure requires updating all JSONB documents
- **No automatic migrations**: Can't use standard migration tools easily

**Example Challenge**:
```json
// Old format
{"control_time_step": 1000}

// New format (renamed field)
{"time_step_ms": 1000}

// Application must handle both
```

### 5. **Debugging and Tooling** ⚠️
- **Harder to inspect**: `psql` queries require JSONB operators
- **Less intuitive**: Developers need to know JSONB syntax
- **Tool limitations**: Some database tools don't display JSONB nicely
- **Query complexity**: Simple questions require complex JSONB queries

**Example**:
```sql
-- Simple question: "What's the average control_time_step?"
-- Requires JSONB extraction
SELECT AVG((input_parameters->>'control_time_step')::int) 
FROM simulations;
```

### 6. **Indexing Overhead** ⚠️
- **GIN indexes are large**: JSONB GIN indexes can be 2-3x larger than B-tree
- **Slower index updates**: GIN indexes are slower to maintain
- **Limited index types**: Can't use specialized indexes (e.g., GiST for ranges)

---

## Alternative: Normalized Columns

### What It Would Look Like
```sql
CREATE TABLE simulations (
    id UUID PRIMARY KEY,
    tle_line1 TEXT NOT NULL,
    tle_line2 TEXT NOT NULL,
    start_date_time TIMESTAMP NOT NULL,
    end_date_time TIMESTAMP NOT NULL,
    control_time_step INTEGER NOT NULL CHECK (control_time_step > 0),
    ground_station_lat DOUBLE PRECISION NOT NULL,
    ground_station_lon DOUBLE PRECISION NOT NULL,
    ground_station_alt DOUBLE PRECISION NOT NULL,
    ground_station_elevation DOUBLE PRECISION NOT NULL,
    status VARCHAR(50),
    -- ... other columns
);
```

### Pros of Normalized Approach
- ✅ **Fast queries**: Direct column indexes, no JSONB overhead
- ✅ **Database constraints**: NOT NULL, CHECK, foreign keys
- ✅ **Type safety**: Native PostgreSQL types
- ✅ **Better tooling**: Standard SQL tools work well
- ✅ **Easier debugging**: Simple column queries

### Cons of Normalized Approach
- ❌ **Schema rigidity**: Adding parameters requires ALTER TABLE
- ❌ **Migration overhead**: Schema changes need careful planning
- ❌ **More columns**: Wider table (currently 9+ columns vs 1 JSONB)
- ❌ **NULL handling**: Optional fields require nullable columns or default values
- ❌ **Atomic updates**: Need to update multiple columns in transaction

---

## Hybrid Approach: JSONB + Generated Columns

### Best of Both Worlds
```sql
CREATE TABLE simulations (
    id UUID PRIMARY KEY,
    input_parameters JSONB NOT NULL,
    
    -- Generated columns for frequently queried fields
    start_date_time TIMESTAMP 
        GENERATED ALWAYS AS (
            (input_parameters->>'start_date_time')::timestamp
        ) STORED,
    control_time_step INTEGER
        GENERATED ALWAYS AS (
            (input_parameters->>'control_time_step')::integer
        ) STORED,
    
    -- Regular columns for metadata
    status VARCHAR(50),
    created_at TIMESTAMP DEFAULT NOW()
);

-- Create indexes on generated columns
CREATE INDEX idx_simulations_start_date 
    ON simulations(start_date_time);
```

### Benefits
- ✅ **Keep JSONB flexibility** for parameter storage
- ✅ **Fast queries** on generated columns with indexes
- ✅ **Database constraints** can be added to generated columns
- ✅ **Best of both worlds**: Flexibility + Performance

### When to Use
- Use when you need to **query/filter by specific parameters**
- Use when you need **database-level constraints** on specific fields
- Use when you need **better tooling support** for specific fields

---

## Recommendation Matrix

| Scenario | Recommendation | Rationale |
|----------|---------------|-----------|
| **Current state** (query only by ID) | ✅ **Keep JSONB** | No performance issues, maximum flexibility |
| **Need to filter by dates** | ⚠️ **Add generated columns** | Enable fast date range queries |
| **Need to filter by ground station** | ⚠️ **Add generated columns** | Enable spatial queries |
| **Need database constraints** | ⚠️ **Add CHECK constraints on JSONB** or **generated columns** | Validate data at DB level |
| **Adding many new parameters** | ✅ **Keep JSONB** | Avoid schema churn |
| **Need to aggregate by parameters** | ⚠️ **Add generated columns** | Enable GROUP BY queries |
| **API returns parameters directly** | ✅ **Keep JSONB** | No transformation needed |

---

## Specific Recommendations for Your Application

### 1. **Keep JSONB for Now** ✅
Your current usage pattern (primary key lookups only) makes JSONB ideal:
- No performance concerns
- Maximum flexibility for future parameters
- Simple application code

### 2. **Add Validation Layer** ⚠️
Enhance `parseInputParameters()` with stricter validation:
```cpp
// Add validation
if (params.control_time_step <= 0) {
    throw std::runtime_error("control_time_step must be positive");
}
if (params.start_date_time >= params.end_date_time) {
    throw std::runtime_error("start_date_time must be before end_date_time");
}
```

### 3. **Consider Generated Columns (Future)** 🔮
If you later need queries like:
- "Find all simulations starting in 2023"
- "Find simulations with ground station in New York"
- "Average control_time_step across all simulations"

Then add generated columns for those specific fields.

### 4. **Add JSONB Index (If Needed)** 🔮
If you need to search within JSONB:
```sql
-- GIN index for JSONB queries
CREATE INDEX idx_simulations_params_gin 
    ON simulations USING GIN (input_parameters);

-- Enables fast queries like:
-- WHERE input_parameters @> '{"control_time_step": 1000}'
```

### 5. **Schema Versioning** 💡
Consider adding a version field:
```json
{
  "schema_version": "1.0",
  "tle_line1": "...",
  ...
}
```
This helps handle schema evolution gracefully.

---

## Performance Comparison

### Current Query Pattern (Primary Key Lookup)
```sql
SELECT input_parameters FROM simulations WHERE id = $1;
```
- **JSONB**: ⚡ **Fast** (indexed primary key, minimal overhead)
- **Normalized**: ⚡ **Fast** (same performance)
- **Verdict**: No difference

### Hypothetical Future Query (Filter by Date)
```sql
-- JSONB approach
SELECT * FROM simulations 
WHERE (input_parameters->>'start_date_time')::timestamp > '2023-01-01';

-- Normalized approach
SELECT * FROM simulations 
WHERE start_date_time > '2023-01-01';
```
- **JSONB**: 🐌 **Slower** (requires JSONB extraction, no index)
- **Normalized**: ⚡ **Fast** (direct column index)
- **With Generated Column**: ⚡ **Fast** (indexed generated column)
- **Verdict**: Normalized or generated column wins

---

## Conclusion

**For your current application: JSONB is the right choice.**

### Keep JSONB Because:
1. ✅ You only query by primary key (no performance penalty)
2. ✅ Parameters may evolve (flexibility is valuable)
3. ✅ Application code is simpler (single JSON structure)
4. ✅ No current need for filtering/aggregation by parameters

### Consider Changes If:
1. ⚠️ You need to filter/query by parameter values → Add generated columns
2. ⚠️ You need database-level validation → Add CHECK constraints or generated columns
3. ⚠️ You need to aggregate by parameters → Add generated columns
4. ⚠️ Performance becomes an issue → Profile and optimize (likely generated columns)

### Migration Path:
If you later need normalized columns, you can:
1. Add generated columns (no data migration needed)
2. Create indexes on generated columns
3. Gradually migrate queries to use generated columns
4. Keep JSONB as source of truth

**Bottom Line**: Your current JSONB approach is appropriate. Don't change it unless you have specific requirements that JSONB can't meet efficiently.

