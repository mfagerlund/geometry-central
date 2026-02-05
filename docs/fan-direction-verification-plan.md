# Fan Direction Verification Plan

## Problem Statement

When VDG (VeryDiscreteGeodesic) takes a deep jump (L2+), the face strip builder must connect the current face to the destination face by "fanning" around a shared vertex. At each such vertex, there are two possible directions:

- **Clockwise fan**: Walk around the vertex one way
- **Counterclockwise fan**: Walk around the vertex the other way

Picking the wrong fan direction creates a sleeve with more faces than necessary. While the funnel algorithm will still find a valid path, the sleeve starts suboptimal and requires more iterations to straighten.

## Goal

Verify that the current implementation picks the optimal fan direction, or quantify how often it picks suboptimally and fix it.

## Verification Approach

### Phase 1: Instrumentation

Add diagnostic code to track fan direction choices:

1. **Locate fan construction code** in `very_discrete_geodesic.cpp`:
   - Find where deep jumps connect faces
   - Identify the decision point for left vs right fan

2. **Add dual-path comparison**:
   - For each fan decision, compute BOTH directions
   - Record: face count, initial funnel distance for each
   - Track which direction was chosen vs which was optimal

### Phase 2: Test Implementation

Create `benchmark/fan_direction_test.cpp`:

```cpp
// For each path:
// 1. Run VDG to get vertex sequence with deep jumps
// 2. At each fan decision point:
//    a. Build sleeve with LEFT fan
//    b. Build sleeve with RIGHT fan
//    c. Compare: face count, initial funnel distance
// 3. Record: correct choices, wrong choices, impact

struct FanDecision {
    size_t pathIndex;
    size_t jumpIndex;
    Vertex fanVertex;

    size_t leftFaceCount;
    size_t rightFaceCount;
    double leftFunnelDist;
    double rightFunnelDist;

    bool choseLeft;
    bool optimalWasLeft;  // based on funnel distance
};

struct FanVerificationResult {
    size_t totalDecisions;
    size_t correctDecisions;
    size_t wrongDecisions;

    double avgExtraFaces;      // when wrong
    double avgExtraDistance;   // when wrong (% of optimal)
    double avgExtraIterations; // when wrong
};
```

### Phase 3: Metrics to Collect

For each mesh (1000 paths):

| Metric | Description |
|--------|-------------|
| Total fan decisions | How many left/right choices were made |
| Correct choices | Chose the direction with shorter funnel distance |
| Wrong choices | Chose the longer direction |
| Wrong choice rate | % of decisions that were suboptimal |
| Avg extra faces | Mean additional faces when wrong |
| Avg extra distance | Mean % longer initial funnel when wrong |
| Avg extra iterations | Mean additional straightening iterations when wrong |

### Phase 4: Test Execution

Run on standard benchmark meshes:

```bash
./fan_direction_test --paths 1000 \
    --mesh bunny-small.obj \
    --mesh spot.obj \
    --mesh stanford-bunny.obj \
    --mesh armadillo.obj \
    --mesh max-planck.obj \
    --mesh happy-buddha.obj
```

Expected output:
```
Fan Direction Verification Results
==================================

Mesh            | Decisions | Correct | Wrong | Wrong% | Extra Faces | Extra Dist
----------------|-----------|---------|-------|--------|-------------|------------
bunny-small     | 2,341     | 2,298   | 43    | 1.8%   | +2.3        | +0.4%
spot            | 3,102     | 3,045   | 57    | 1.8%   | +1.9        | +0.3%
...

Overall: X% wrong choices, averaging Y extra faces and Z% longer initial distance
```

## Decision Tree

Based on results:

### If wrong choice rate < 1%:
- Document in paper: "Fan direction is selected optimally in >99% of cases"
- No code changes needed

### If wrong choice rate 1-5%:
- Investigate the failure cases
- Determine if there's a pattern (certain mesh types, path lengths, etc.)
- Consider if fixing is worth the complexity

### If wrong choice rate > 5%:
- This is a bug - the heuristic is broken
- Fix the fan direction selection logic
- Re-run benchmarks to show improvement

## Potential Fixes (if needed)

### Option A: Try-both approach
At each fan decision, compute both directions and pick the better one.
- Pro: Always optimal
- Con: 2x cost at fan points (but fan building is cheap vs funnel)

### Option B: Geometric heuristic
Use the target vertex position to inform direction choice:
- Compute angle from current position to target
- Pick fan direction that "opens toward" the target
- Pro: O(1) decision
- Con: May not always be optimal

### Option C: Distance-based heuristic
Pick the fan direction where the first face is closer to the target:
- Pro: Simple, often correct
- Con: Greedy, may fail on curved paths

## Paper Addition (if verification passes)

Add to "Sleeve Invariant Verification" section:

### Fan Direction Optimality

When building face strips from VDG deep jumps, the algorithm must choose a fan direction around shared vertices. We verified this choice is optimal:

| Mesh | Fan Decisions | Optimal | Suboptimal | Rate |
|------|---------------|---------|------------|------|
| ... results ... |

The fan direction heuristic achieves >X% optimality. Suboptimal choices add an average of Y faces to the sleeve, resulting in Z% longer initial funnel distance and W additional straightening iterations.

## Implementation Checklist

- [x] Read `very_discrete_geodesic.cpp` to find fan construction logic
- [x] Identify the left/right decision point
- [x] Create `fan_direction_test.cpp` with dual-path comparison
- [x] Add to `benchmark/CMakeLists.txt`
- [x] Run on 6 benchmark meshes (1000 paths each)
- [x] Analyze results
- [x] If >5% wrong: fix the heuristic
- [ ] Update paper with verification results
- [ ] Commit and push

## Actual Results

### Before Fix (geometric heuristic in `handleEdgeStep`)

| Mesh | Decisions | Optimal | Suboptimal | Rate | Avg Extra Angle |
|------|-----------|---------|------------|------|-----------------|
| bunny-small | 216 | 116 | 100 | 53.7% | +29.8 deg |
| spot | 1,071 | 835 | 236 | 78.0% | +60.0 deg |
| stanford-bunny | 221 | 99 | 122 | 44.8% | +25.5 deg |
| armadillo | 2,195 | 1,393 | 802 | 63.5% | +47.7 deg |
| max-planck | 1,703 | 1,182 | 521 | 69.4% | +48.5 deg |
| happy-buddha | 369 | 213 | 156 | 57.7% | +38.1 deg |

**Overall: 66.5% optimal**, 33.5% suboptimal with +46.3 deg extra angle on average.

### After Fix (try both directions, pick smaller angle sum)

| Mesh | Decisions | Optimal | Suboptimal | Rate |
|------|-----------|---------|------------|------|
| bunny-small | 216 | 216 | 0 | 100.0% |
| spot | 1,071 | 1,071 | 0 | 100.0% |
| stanford-bunny | 221 | 221 | 0 | 100.0% |
| armadillo | 2,195 | 2,195 | 0 | 100.0% |
| max-planck | 1,703 | 1,703 | 0 | 100.0% |
| happy-buddha | 369 | 369 | 0 | 100.0% |

**Overall: 100% optimal** across 5,775 fan decisions.

### Fix Applied

Modified `handleEdgeStep()` in `very_discrete_geodesic.cpp` to:
1. Try both CW and CCW directions
2. Compute angle sum for each walk
3. Pick the direction with smaller angle sum

This matches the approach already used in `handleApexJump()`.

## Files to Modify/Create

1. **Create**: `benchmark/fan_direction_test.cpp`
2. **Modify**: `benchmark/CMakeLists.txt` (add new target)
3. **Possibly modify**: `src/surface/very_discrete_geodesic.cpp` (if fix needed)
4. **Update**: `paper-draft.md` (add verification results)
