# CLAUDE.md - geometry-central GFR Fork

## Build Commands

CMake is at: `"/c/Program Files/CMake/bin/cmake.exe"` (or add to PATH)

```bash
# Configure (from geometry-central-gfr root)
cd build
"/c/Program Files/CMake/bin/cmake.exe" ..

# Build
"/c/Program Files/CMake/bin/cmake.exe" --build . --config Release

# Or use the alias
CMAKE="/c/Program Files/CMake/bin/cmake.exe"
```

## Test Commands

```bash
# Run tests (after building test project)
cd test/build
"$CMAKE" ..
"$CMAKE" --build . --config Release
./bin/Release/geometry-central-test.exe
```

## GFR Implementation

Funnel Geodesics (GFR) - alternative to FlipOut for computing geodesic paths.

**Files:**
- `include/geometrycentral/surface/funnel_geodesics.h` - Header
- `src/surface/funnel_geodesics.cpp` - Implementation
- `test/src/funnel_geodesics_test.cpp` - Tests
- `benchmark/geodesic_comparison.cpp` - GFR vs FlipOut comparison

**API:**
```cpp
auto path = computeFunnelGeodesic(mesh, geometry, startVert, endVert);
double len = path->length();
auto points = path->getPath();  // std::vector<SurfacePoint>
```

**Algorithm:** FaceStripBuilding (Dijkstra corridor) → Funnel (unfold to 2D) → IterativeStraightener (quality-gated corner flipping)

## Current Status

**Working:**
- Phase 1: 2D Flattening (flattenSleeve)
- Phase 2: Portal Building (buildPortals)
- Phase 3: Funnel Algorithm (runFunnel)
- Phase 5: Corner Analysis (analyzeCorners)
- Phase 6: Flip Action Computation (computeFlipAction)
- Phase 7: Apply Flip (applyFlip)
- Phase 8: Iterative Straightening loop

**Needs Work:**
- Phase 4: buildFaceStrip has connectivity issues
  - Some face strips have gaps (non-adjacent consecutive faces)
  - Need to port walk-based approach from C# FaceStripWalker.cs
  - Current implementation uses simple edge-face collection

## Submodules

If build fails with missing happly.h:
```bash
git submodule update --init --recursive
```
