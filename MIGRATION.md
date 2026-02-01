# GFR Port Migration Notes

## Overview

This fork adds Funnel Geodesics (GFR) alongside FlipOut for comparison.

## TODO: Migrate from Colonel

### 1. Comparison Data Generation

Port the batch comparison logic from:
- `C:\Dev\Colonel\Colonel.Meshing\GreedyFunnelRefinement\FlipOutComparisonTests.cs`

This generates vertex pairs and runs both algorithms for comparison.

### 2. Test Meshes

Copy test meshes to `test/assets/` or reference from a shared location:
- Bunny, Spot, Pig (small)
- StanfordBunny, Armadillo, MaxPlanck, HappyBuddha (large)

Current location: Check Colonel test assets or download from common sources.

### 3. FlipOut Batch Runner

The existing `C:\Dev\flip-geodesics-demo\src\flipout-batch.cpp` has batch testing code.
Migrate the relevant parts here as a standalone benchmark tool.

## File Structure

```
geometry-central-gfr/
├── include/geometrycentral/surface/
│   ├── flip_geodesics.h          # Existing FlipOut
│   └── funnel_geodesics.h        # NEW: GFR header
├── src/surface/
│   ├── flip_geodesics.cpp        # Existing FlipOut
│   └── funnel_geodesics.cpp      # NEW: GFR implementation
├── test/src/
│   └── funnel_geodesics_test.cpp # NEW: GFR tests
└── benchmark/
    └── geodesic_comparison.cpp   # NEW: FlipOut vs GFR benchmark
```

## Port Plan

See `docs/GFR_PORT_PLAN.md` for the phased implementation plan.

## Branch Strategy

- `funnel-geodesics` - main development branch
- Keep in sync with upstream via `git fetch upstream && git merge upstream/master`
