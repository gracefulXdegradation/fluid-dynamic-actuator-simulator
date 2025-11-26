# Database Access Approach: Trade-off Analysis

## Executive Summary

**Recommendation: Continue with Direct SQL (current approach) with minor improvements**

For this orbital mechanics simulation application, **direct SQL queries** are the optimal choice. The data model is simple, queries are straightforward, and the application requires PostgreSQL-specific features and high-performance batch operations that ORMs struggle with.

---

## Application Context

### Data Model Characteristics
- **Schema Complexity**: Very simple - only 2 tables
  - `simulations`: Metadata table (id, status, input_parameters JSONB, timestamps)
  - `simulation_metrics`: Time-series data (simulation_id FK, step_index, timestamp, arrays of vectors/matrices)
- **Relationships**: Minimal (single foreign key relationship)
- **Query Patterns**: Simple CRUD operations
  - SELECT by ID
  - UPDATE by ID
  - Batch INSERT with ON CONFLICT
  - No complex joins, aggregations, or subqueries

### Application Characteristics
- **Language**: C++ (performance-critical)
- **Domain**: Scientific computing (orbital mechanics simulation)
- **Data Volume**: High-volume time-series writes (thousands of metrics per simulation)
- **Database**: PostgreSQL with heavy use of:
  - JSONB for flexible parameter storage
  - Array types for Eigen vectors/matrices
  - PostgreSQL-specific syntax (ON CONFLICT, array literals)

### Current Implementation
- Uses `pqxx` (PostgreSQL C++ library)
- Direct SQL with parameterized queries
- Manual conversion of Eigen types to PostgreSQL array format
- Batch inserts in transactions

---

## Option 1: ORM (Object-Relational Mapping)

### Pros
- ✅ Type safety and compile-time checking
- ✅ Automatic object-relational mapping
- ✅ Reduces boilerplate for simple CRUD
- ✅ Database-agnostic (in theory)

### Cons
- ❌ **Poor C++ ORM ecosystem**: Limited mature options (ODB, QxOrm, SOCI) - none are as mature as Python/Ruby ORMs
- ❌ **Performance overhead**: ORM abstraction layers add overhead for high-volume operations
- ❌ **PostgreSQL-specific features**: Most C++ ORMs struggle with:
  - JSONB operations
  - Array types (especially custom array formats)
  - ON CONFLICT syntax
  - Custom type conversions (Eigen → PostgreSQL arrays)
- ❌ **Learning curve**: Team needs to learn ORM-specific syntax and patterns
- ❌ **Limited flexibility**: Hard to optimize batch inserts or use advanced PostgreSQL features
- ❌ **Eigen integration**: Would require custom type adapters for every Eigen type
- ❌ **Overkill for simple schema**: Your 2-table schema doesn't benefit from ORM complexity

### Example Issues You'd Face
```cpp
// With ORM, you'd need custom adapters for:
Eigen::Vector3d → PostgreSQL array[3]
Eigen::Vector4d → PostgreSQL array[4]
Eigen::Matrix<double, 15, 1> → PostgreSQL array[15]
JSONB operations → Custom serialization

// Your current approach is cleaner:
std::string euler_angles_str = "{" + 
    std::to_string(metric.euler_angles(0)) + "," +
    std::to_string(metric.euler_angles(1)) + "," +
    std::to_string(metric.euler_angles(2)) + "}";
```

### Verdict: ❌ **Not Recommended**
The C++ ORM ecosystem is weak, and your use case requires PostgreSQL-specific features that ORMs handle poorly.

---

## Option 2: Query Builder

### Pros
- ✅ Type-safe query construction
- ✅ Compile-time query validation (some libraries)
- ✅ More readable than raw SQL strings
- ✅ Protection against SQL injection
- ✅ Easier to compose dynamic queries

### Cons
- ❌ **Limited C++ options**: Few mature query builders for C++/PostgreSQL
  - Most are C++ wrappers around SQL strings anyway
- ❌ **PostgreSQL-specific features**: Still need to drop to raw SQL for:
  - JSONB operations
  - Array type handling
  - ON CONFLICT syntax
- ❌ **Eigen integration**: Still requires manual conversion
- ❌ **Additional dependency**: Another library to maintain
- ❌ **Learning curve**: Team needs to learn query builder API
- ❌ **Performance**: Adds abstraction layer overhead

### Example with Query Builder
```cpp
// Hypothetical query builder (most C++ ones don't support this well)
auto query = db.select("input_parameters")
    .from("simulations")
    .where("id", "=", simulation_id);

// But you'd still need raw SQL for:
// - JSONB operations
// - Array inserts
// - ON CONFLICT
```

### Verdict: ⚠️ **Marginally Useful**
Query builders could help with basic queries, but you'd still need raw SQL for PostgreSQL-specific features. The benefit is minimal for your simple query patterns.

---

## Option 3: Direct SQL (Current Approach)

### Pros
- ✅ **Full PostgreSQL feature access**: JSONB, arrays, ON CONFLICT, etc.
- ✅ **Performance**: No abstraction overhead, direct control
- ✅ **Eigen integration**: Straightforward manual conversion (already working)
- ✅ **Mature library**: `pqxx` is well-established and reliable
- ✅ **Simple queries**: Your queries are straightforward, SQL is readable
- ✅ **Batch operations**: Full control over transaction batching
- ✅ **Debugging**: Easy to copy SQL to `psql` for testing
- ✅ **No additional dependencies**: Already using `pqxx`

### Cons
- ❌ **String concatenation**: Manual SQL string building (mitigated by parameterized queries)
- ❌ **No compile-time validation**: SQL errors only at runtime
- ❌ **Boilerplate**: Some repetitive code for array formatting
- ❌ **SQL injection risk**: Mitigated by parameterized queries (`$1`, `$2`, etc.)

### Current Implementation Quality
Your current code already uses best practices:
- ✅ Parameterized queries (`$1`, `$2`) - prevents SQL injection
- ✅ Transactions for batch operations
- ✅ Error handling
- ✅ Connection management

### Verdict: ✅ **Recommended - Continue with Improvements**

---

## Recommended Approach: Enhanced Direct SQL

### Keep Current Approach, Add These Improvements:

#### 1. **Helper Functions for Array Formatting**
```cpp
// In Database.cpp or a helper header
namespace db_helpers {
    template<typename Derived>
    std::string eigenVectorToArray(const Eigen::MatrixBase<Derived>& vec) {
        std::ostringstream oss;
        oss << "{";
        for (int i = 0; i < vec.size(); ++i) {
            if (i > 0) oss << ",";
            oss << vec(i);
        }
        oss << "}";
        return oss.str();
    }
}

// Usage:
std::string euler_angles_str = db_helpers::eigenVectorToArray(metric.euler_angles);
```

#### 2. **Prepared Statements for Repeated Queries**
```cpp
// In Database class initialization
void Database::prepareStatements() {
    conn_->prepare("get_simulation_params",
        "SELECT input_parameters FROM simulations WHERE id = $1");
    conn_->prepare("update_status_running",
        "UPDATE simulations SET status = $1, started_at = CURRENT_TIMESTAMP WHERE id = $2");
    // ... etc
}

// Usage:
pqxx::result result = txn.exec_prepared("get_simulation_params", simulation_id);
```

#### 3. **Batch Insert Optimization**
Consider using `COPY` for very large batches:
```cpp
// For 10,000+ rows, COPY is faster than INSERT
void Database::writeMetricsCopy(const std::string& simulation_id, 
                                 const std::vector<SimulationMetric>& metrics) {
    // Use PostgreSQL COPY command for bulk inserts
    // Requires formatting data as CSV or binary
}
```

#### 4. **Type-Safe Query Wrapper (Optional)**
Create a thin wrapper for common patterns:
```cpp
template<typename T>
T Database::queryOne(const std::string& sql, const std::vector<std::string>& params) {
    // Type-safe single-row query wrapper
}
```

---

## Comparison Matrix

| Criteria | ORM | Query Builder | Direct SQL (Enhanced) |
|----------|-----|---------------|----------------------|
| **C++ Ecosystem Maturity** | ⭐⭐ Poor | ⭐⭐⭐ Limited | ⭐⭐⭐⭐⭐ Excellent (pqxx) |
| **PostgreSQL Features** | ⭐⭐ Limited | ⭐⭐⭐ Partial | ⭐⭐⭐⭐⭐ Full access |
| **Performance** | ⭐⭐⭐ Good | ⭐⭐⭐⭐ Very Good | ⭐⭐⭐⭐⭐ Excellent |
| **Eigen Integration** | ⭐⭐ Difficult | ⭐⭐⭐ Moderate | ⭐⭐⭐⭐ Easy |
| **Code Readability** | ⭐⭐⭐⭐ Good | ⭐⭐⭐⭐ Good | ⭐⭐⭐ Moderate |
| **Learning Curve** | ⭐⭐ Steep | ⭐⭐⭐ Moderate | ⭐⭐⭐⭐ Easy |
| **Maintenance** | ⭐⭐ Complex | ⭐⭐⭐ Moderate | ⭐⭐⭐⭐ Simple |
| **Type Safety** | ⭐⭐⭐⭐⭐ Excellent | ⭐⭐⭐⭐ Very Good | ⭐⭐⭐ Good |
| **Flexibility** | ⭐⭐ Limited | ⭐⭐⭐ Moderate | ⭐⭐⭐⭐⭐ Full |
| **Your Use Case Fit** | ⭐⭐ Poor | ⭐⭐⭐ Moderate | ⭐⭐⭐⭐⭐ Excellent |

---

## Decision Factors Summary

### Why Direct SQL Wins for Your Application:

1. **Simple Schema**: 2 tables don't justify ORM complexity
2. **PostgreSQL-Specific Features**: Heavy use of JSONB, arrays, ON CONFLICT
3. **Performance-Critical**: High-volume batch inserts need direct control
4. **Eigen Integration**: Manual conversion is straightforward and explicit
5. **Mature Library**: `pqxx` is battle-tested and well-documented
6. **Team Familiarity**: SQL is more universally understood than ORM APIs
7. **Debugging**: Easy to test queries in `psql`

### When to Reconsider:

- **If schema grows complex**: 10+ tables with many relationships → consider ORM
- **If queries become complex**: Many joins, aggregations, subqueries → consider query builder
- **If team grows**: More developers unfamiliar with SQL → consider query builder for safety
- **If database changes**: Need to support multiple databases → consider ORM/query builder

---

## Conclusion

**Stick with Direct SQL**, but enhance it with:
1. Helper functions for array formatting (reduce boilerplate)
2. Prepared statements (improve performance for repeated queries)
3. Consider COPY for very large batches (10,000+ rows)

Your current approach is appropriate for this application. The simplicity of your data model and the need for PostgreSQL-specific features make direct SQL the pragmatic choice. Focus on small improvements rather than architectural changes.

---

## References

- [pqxx Documentation](https://libpqxx.readthedocs.io/)
- [PostgreSQL Arrays](https://www.postgresql.org/docs/current/arrays.html)
- [PostgreSQL JSONB](https://www.postgresql.org/docs/current/datatype-json.html)
- [C++ ORM Comparison](https://github.com/topics/cpp-orm)

