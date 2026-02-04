# Portal Passability Discrepancy Investigation

## The Problem

GFR paths became LONGER (worse) after changes were made to align VDG's portal checking with flattenSleeve's. This suggests **the wrong component was modified**. The original VDG deep jumps were finding genuinely shorter paths.

**Before "fix":** VDG found deep jumps → shorter paths than FlipOut
**After "fix":** VDG and flattenSleeve agree → paths are LONGER than FlipOut

## The Core Issue

VeryDiscreteGeodesic (VDG/L5) explores from a vertex and determines which deep vertices (L2, L3, L4, L5) are reachable via straight lines through portals. It uses 2D unfolding to check portal crossings.

When the face strip is built and flattened by `flattenSleeve`, the funnel algorithm should find the same straight-line path. But sometimes **VDG says a path is valid, and flattenSleeve/funnel disagrees**.

This is a fundamental correctness bug. Both are doing isometric unfolding of triangles. If done correctly, they MUST agree on whether a straight line crosses a portal.

## What We Need To Find

**A concrete case where:**
1. VDG's `explore()` says vertex X is reachable from corner C (straight line crosses all portals)
2. The face strip is built containing faces [f0, f1, f2, ...]
3. `flattenSleeve` flattens those faces
4. The funnel algorithm (or manual inspection) shows the straight line does NOT cross the portals

## Key Questions

1. Are VDG and flattenSleeve unfolding the **same faces**?
2. Are they unfolding faces in the **same order**?
3. Are the **portal edges** the same in both?
4. Are the **2D positions** of portal endpoints the same (or at least producing same sidedness)?

## Coordinate Systems (for reference)

**VDG explore():**
```cpp
// Portal centered on X-axis
Vector2 flat_heL1_a = {-portalLen / 2.0, 0.0};
Vector2 flat_heL1_b = {portalLen / 2.0, 0.0};
// v0 on positive Y, targets on negative Y
Vector2 flat_v0 = computeTriangleApex(..., pickPositiveY: true);
Vector2 flat_L1 = computeTriangleApex(..., pickPositiveY: false);
```

**C# FlattenFaceStripD:**
```csharp
// Entry at origin, edge along +X
var p0 = Vec2D.Zero;
var p1 = new Vec2D(d01, 0);
// New vertices placed opposite "reference point" (fell-off vertex)
var newPos = CircleIntersectD(h1Pos, r1, h2Pos, r2, d, refPos);
```

**Both approaches are valid** for isometric unfolding. The "fell off" reference point approach is correct - always unfold outward, new vertex opposite the reference.

## Reproducing on a Flat Mesh

**THIS WOULD HELP IMMENSELY.**

On a flat mesh (e.g., a planar triangulated grid):
- All faces lie in a plane
- Unfolding should produce the EXACT same positions as 3D (no distortion)
- Any discrepancy between VDG and flattenSleeve would be clearly a bug, not numerical issues

**Test approach:**
1. Create/use a flat triangulated mesh (e.g., `FlatSquare` or similar)
2. Pick a vertex pair where VDG uses a deep jump (L3+)
3. Log VDG's portal positions and crossing checks
4. Log flattenSleeve's portal positions
5. Compare - find where they disagree

## Files Involved

- `src/surface/very_discrete_geodesic.cpp` - VDG explore(), portal crossing checks
- `src/surface/funnel_geodesics.cpp` - flattenSleeve(), buildPortals()
- `include/geometrycentral/surface/very_discrete_geodesic.h` - data structures

## C# Reference (working implementation)

- `Colonel.Meshing/GreedyFunnelRefinement/DiscreteGeodesics/VeryDiscreteGeodesicExplorer.cs`
- `Colonel.Meshing/Funnels/FunnelAlgorithm.cs` - FlattenFaceStripD()

## Next Steps

1. **Create diagnostic tool** that logs both VDG's and flattenSleeve's portal computations for a single path
2. **Test on flat mesh** to eliminate numerical issues
3. **Find first discrepancy** - which portal, which face, what positions
4. **Determine root cause** - wrong faces? wrong order? wrong edge? wrong side calculation?
5. **Fix the correct component** - likely flattenSleeve needs to match VDG, not vice versa

## Important Note

The uncommitted changes to `flattenSleeve` (making it use "centered portal like L5") may need to be reverted. The original C# convention was working. The fix should make the components agree, but we changed the wrong one.

```bash
# To see uncommitted changes:
git diff HEAD -- src/surface/funnel_geodesics.cpp

# To revert if needed:
git checkout -- src/surface/funnel_geodesics.cpp
```
