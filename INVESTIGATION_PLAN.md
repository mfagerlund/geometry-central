# C++ GFR Investigation Plan

## Current Status (2026-02-01) - FULLY RESOLVED + CLEANED UP

### Original Issue
The C++ GFR port had a -55% outlier that was **physically impossible** (shorter than the exact geodesic).

### Root Cause
The `convertPath` function in `very_discrete_geodesic.cpp` was a simplified port that missed critical logic from the C# `FaceStripWalker.ConvertPath`:
1. Missing `RemoveRange` logic for backtracking (truncating face list when revisiting a face)
2. Missing `CountBadFaces` to prefer walks avoiding entry/exit vertices
3. Missing `finalIdx` check in edge step handling

This caused disconnected face strips, which led to invalid flattening and impossible path lengths.

### Fixes Applied
1. Rewrote `convertPath` as a proper line-by-line port of `FaceStripWalker.ConvertPath`
2. Added `handleApexJump` and `handleEdgeStep` helper functions matching C#
3. Added connectivity validation as a fallback (falls back to Dijkstra if strip is disconnected)
4. Implemented proper phase-based iterative straightening matching C# `InterativeStraightener.StraightenWithPhases()`
5. **Added `CachedVeryDiscreteGeodesicPathfinder`** with aggressive caching (see Performance section)

### Cleanup (2026-02-01)
Removed all investigation/research code to keep only the core algorithm:
- Fixed exploration depth to **L5 only** (11 candidates: L1, L2L-L5L, L2R-L5R, L5LM, L5RM)
- Removed configurable `ExplorationDepth` enum and `VeryDiscreteGeodesicConfig`
- Removed L6/L7 candidate handling (enum values, computation, face-crossing functions)
- Removed `enableDirectionPruning` and `directionThresholds` (incompatible with caching)
- Single straightening strategy (corner greedy) - no multi-strategy investigation code
- Hardcoded `HEURISTIC_WEIGHT = 1.0` (standard A*)

### Final Results

| Metric | C++ GFR | FlipOut | C# GFR | Notes |
|--------|---------|---------|--------|-------|
| Time | **~1221 ms** | 3476 ms | 1848 ms | **C++ 2.85x faster than FlipOut** |
| Mean diff | -0.10% | baseline | -0.12% | GFR finds slightly shorter paths |
| Total dist | 2602.79 | 2605.77 | 2602.67 | GFR wins by 0.11% |
| Min diff | -7.76% | - | - | **Legitimate** - GFR closer to exact |
| GFR shorter | 202 (20.2%) | - | - | |
| FlipOut shorter | 124 (12.4%) | - | - | |
| Match (<0.01%) | 674 (67.4%) | - | - | |

All paths are now valid (>= MMP exact distance).

#### Timing Breakdown (1000 Bunny paths, L5 exploration)

| Phase | C++ Time | Per-path |
|-------|----------|----------|
| A* pathfinding | ~1045 ms | 1.05 ms |
| Flatten/Funnel | ~10 ms | 0.01 ms |
| Straightening | ~164 ms | 0.16 ms |
| **Total** | **~1221 ms** | **1.22 ms** |

#### Optimizations Applied (2026-02-01)

1. **flattenSleeve O(n²) → O(n)**: Replaced O(n²) vertex lookup scan with `VertexData<char>` known mask
2. **Phase-stamped rejection array**: Replaced `std::set<size_t>` with `std::vector<uint32_t>` phase stamps for O(1) lookups

Combined effect: 44% reduction in total time (2180ms → 1221ms)

### C++ vs C# GFR Comparison

| Metric | Value |
|--------|-------|
| MATCH (<0.001%) | 7/10 (70%) |
| CLOSE (<0.1%) | 2/10 (20%) |
| DIFF (>0.1%) | 1/10 (10%) |
| Mean |diff| | 0.0167% |
| Max |diff| | 0.1031% |

Small differences are expected due to:
- Different initial face strips from VeryDiscreteGeodesic
- Floating point precision (C# float vs C++ double)
- Full recomputation vs partial recomputation in C#

## MMP Validation (Outlier Analysis)

The worst GFR vs FlipOut difference (-7.03%) was validated:

| Algorithm | Distance | vs MMP |
|-----------|----------|--------|
| MMP (exact) | 1.76544 | baseline |
| GFR | 1.77952 | +0.80% |
| FlipOut | 1.91401 | +8.42% |

**Conclusion**: GFR finds a significantly better path (closer to exact geodesic).

## Performance Investigation - RESOLVED

### Original Problem
C++ GFR took 4.17 ms/path vs FlipOut's 3.48 ms/path (17% slower)

### Root Cause
Missing caching that C# has:

1. **Per-corner exploration cache** (biggest impact)
   - C# `CachedGeodesicPathfinder._explorationCache` caches `ExplorationResult` per corner
   - C++ was recomputing L1-L5 unfolding for every corner in every path query

2. **Per-vertex boundary edge count**
   - C# `Vertex.BoundaryEdgeCountCached` computes once per vertex
   - C++ was iterating halfedges for every neighbor check

3. **Per-vertex adjacency arrays**
   - C# `AdjacentVerticesCached`, `AdjacentCornersCached` computed once
   - C++ created new iterators each call

4. **Container pooling**
   - C# uses `[ThreadStatic]` pools for A* containers
   - C++ allocated fresh containers per query

### Fixes Applied

1. **CachedVeryDiscreteGeodesicPathfinder** - Per-corner exploration caching
   - Caches `ExplorationResult` per corner (keyed by corner index)
   - Pre-computes boundary edge counts for all vertices at construction
   - Caches adjacent vertices and corners arrays per vertex
   - Reuses A* containers across queries
   - 98.7% cache hit rate in benchmarks

2. **Return by const reference** - Avoid copying large structs
   - `getOrComputeExploration()` returns `const ExplorationResult&`
   - `getNeighbors()` returns `const std::vector<...>&`
   - PathStep stored by const reference in iteration

3. **Move semantics** - Avoid copying PathStep with crossedFaces
   - `std::move(step)` when pushing to neighbors vector

4. **Default exploration level** - Changed to L5 (matching C# benchmark)

### Performance Results

| Optimization | A* Time | Total Time | Speedup vs FlipOut |
|--------------|---------|------------|--------------------|
| Original (no caching) | ~4000 ms | ~4171 ms | 0.83x (slower) |
| + Exploration caching | ~2557 ms | ~3100 ms | 1.12x |
| + Return by reference | ~2401 ms | ~3000 ms | 1.16x |
| + Neighbors by reference | ~1810 ms | ~2564 ms | 1.36x |
| + Move PathStep | **~1430 ms** | **~2180 ms** | **1.59x** |

C# GFR achieves 1.88x speedup (1848 ms). Remaining gap (~330 ms) is likely from:
- C# `[ThreadStatic]` memory pooling
- Different hash map implementations
- JIT inlining optimizations

## Files Involved

- `src/surface/funnel_geodesics.cpp` - Main GFR implementation + iterative straightener
- `src/surface/very_discrete_geodesic.cpp` - A* pathfinder + FaceStripWalker port + **CachedVeryDiscreteGeodesicPathfinder**
- `include/geometrycentral/surface/very_discrete_geodesic.h` - **CachedVeryDiscreteGeodesicPathfinder** class declaration
- `test/src/funnel_geodesics_test.cpp` - Comparison tests + MMP validation
- `include/geometrycentral/surface/exact_geodesics.h` - Ground truth validation
