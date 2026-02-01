# GFR C++ Port Plan

Bottom-up port of Greedy Funnel Refinement to geometry-central.
Each phase is testable in isolation before proceeding.

## Phase 0: Project Setup
**Goal**: Build environment, file structure, basic test harness

**Tasks**:
- [ ] Create `funnel_geodesics.h` / `funnel_geodesics.cpp` in geometry-central
- [ ] Add to CMakeLists.txt
- [ ] Create test file `test_funnel_geodesics.cpp`
- [ ] Verify builds and links

**Test**: Empty test compiles and runs

---

## Phase 1: 2D Flattening
**Goal**: Unfold a face strip to 2D plane

**Port from**: `FunnelAlgorithm.FlattenFaceStripD()` in Sleeve.cs

**Data structures**:
```cpp
struct FlattenedSleeve {
  std::vector<Face> faces;
  VertexData<Vector2> flatPos;  // 2D positions
  Vertex entry, exit;
};
```

**Algorithm**:
1. Place entry vertex at origin
2. For each face in strip, unfold next vertex using edge lengths + angle
3. Use double precision for accuracy

**Test cases**:
- Single triangle: 3 vertices at known positions
- Two adjacent triangles: verify shared edge preserved
- Long strip (10 faces): verify no accumulated drift vs known mesh

**Validation**: Compare `flatPos` distances to 3D edge lengths (should match exactly)

---

## Phase 2: Portal Construction
**Goal**: Build portal edges between consecutive faces

**Port from**: `Sleeve.BuildPortalsList()`

**Data structure**:
```cpp
struct Portal {
  Vector2 left, right;  // 2D positions of portal endpoints
  Vertex leftVert, rightVert;
};
```

**Algorithm**:
1. For each pair of adjacent faces, find shared edge
2. Determine left/right based on sleeve direction
3. Store 2D positions from Phase 1

**Test cases**:
- Two triangles: one portal, verify left/right orientation
- Strip of N faces: N-1 portals, all properly oriented

---

## Phase 3: Funnel Algorithm
**Goal**: Find shortest path through portal sequence

**Port from**: `FunnelAlgorithm.Populate()`

**Algorithm** (Lee-Preparata):
1. Initialize funnel with entry point as apex
2. For each portal, update left/right funnel edges
3. When funnel inverts, emit waypoint and restart
4. Connect to exit

**Output**:
```cpp
struct FunnelResult {
  std::vector<Vector2> waypoints2D;
  std::vector<Vertex> waypointVertices;  // mesh vertices at corners
  double distance;
};
```

**Test cases**:
- Straight corridor (no turns): distance = Euclidean entry→exit
- Single left turn: one waypoint
- Single right turn: one waypoint
- S-curve: two waypoints
- Already-straight path: 0 interior waypoints

**Validation**: Path distance ≤ discrete edge-walk distance (always)

---

## Phase 4: Face Strip Building
**Goal**: Build initial sleeve from vertex pair using Dijkstra

**Port from**: `FaceStripBuilder.FindSleeve()`

**Uses**: geometry-central's `shortestEdgePath()` (already exists)

**Algorithm**:
1. Run Dijkstra/A* to get vertex path
2. Convert vertex path to face strip (faces touching the path)
3. Deduplicate and order faces

**Output**: `FlattenedSleeve` from Phase 1

**Test cases**:
- Adjacent vertices: 1-2 faces
- Vertices on same face: 1 face
- Opposite sides of mesh: valid connected strip

**Validation**:
- First face contains entry vertex
- Last face contains exit vertex
- Consecutive faces share an edge

---

## Phase 5: Corner Analysis
**Goal**: Identify which path corners can be improved

**Port from**: `WaypointCornerAnalyzer.Populate()`

**Data structure**:
```cpp
struct WaypointCorner {
  Vertex vertex;
  size_t faceIndex;      // position in sleeve
  double angleError;     // degrees from straight (180°)
  bool wantsToFlip;      // angleError > threshold
};
```

**Algorithm**:
1. For each interior waypoint vertex
2. Compute angle at vertex in 2D
3. If angle < 180° - epsilon, corner wants to flip

**Test cases**:
- Straight path: no corners want to flip
- Single corner at 170°: wants to flip (10° error)
- Corner at 179°: doesn't want to flip (below threshold)

---

## Phase 6: Flip Action Computation
**Goal**: Compute which faces to add/remove for a corner flip

**Port from**: `WaypointCornerFlipAction.Compute()`

**Data structure**:
```cpp
struct FlipAction {
  std::vector<Face> removeFaces;
  std::vector<Face> addFaces;
  bool canFlip;  // false if hits mesh boundary
};
```

**Algorithm**:
1. At corner vertex, find the "wedge" of faces on the short side
2. Remove faces currently in sleeve within wedge
3. Add faces on opposite side of wedge
4. Check boundary constraints

**Test cases**:
- Interior vertex: can flip, faces swap correctly
- Boundary vertex: canFlip = false
- Vertex with one adjacent face: canFlip = false

---

## Phase 7: Flip Application
**Goal**: Apply a flip action to modify the sleeve

**Port from**: `Sleeve.ApplyFlip()` (simple rebuild path, not partial)

**Algorithm**:
1. Remove `removeFaces` from sleeve
2. Insert `addFaces` at same position
3. Rebuild flat positions (full recompute)
4. Rebuild portals

**Test cases**:
- Flip a corner, verify face list changed correctly
- Flip and re-run funnel, verify shorter path
- Flip back (reverse action), verify original state restored

---

## Phase 8: Iterative Straightener (Core Loop)
**Goal**: Quality-gated corner flipping until convergence

**Port from**: `InterativeStraightener.StraightenWithPhases()` ONLY

**Algorithm**:
```
rejected = {}
while iteration < maxIters:
    corner = selectMostAcuteCorner(excluding rejected)
    if corner is None: break

    action = computeFlipAction(corner)
    if not action.canFlip:
        rejected.add(corner.vertex)
        continue

    oldDistance = path.distance
    sleeve = applyFlip(action)
    newPath = runFunnel(sleeve)

    if newPath.distance < oldDistance - epsilon:
        # Good flip - keep it
        path = newPath
    else:
        # Bad flip - reverse and reject
        sleeve = applyFlip(action.reverse())
        rejected.add(corner.vertex)
```

**Test cases**:
- Already optimal path: 0 iterations
- Single improvable corner: 1 iteration
- Path requiring multiple phases: correct convergence
- Stress test: 1000 random vertex pairs, all converge

---

## Phase 9: End-to-End API
**Goal**: Clean public API matching FlipOut style

**API**:
```cpp
class FunnelGeodesic {
public:
  // Construct and compute geodesic
  static std::unique_ptr<FunnelGeodesic> compute(
      ManifoldSurfaceMesh& mesh,
      VertexPositionGeometry& geom,
      Vertex start, Vertex end);

  // Query results
  double length() const;
  std::vector<SurfacePoint> getPath() const;
  size_t iterationCount() const;
};
```

**Test cases**:
- Compare to FlipOut on same vertex pairs
- Verify path length ≤ FlipOut path length (or within epsilon)
- Benchmark: 1000 paths, measure time vs FlipOut

---

## Phase 10: Validation & Benchmarking
**Goal**: Comprehensive comparison with FlipOut

**Tests**:
- [ ] Same 7 meshes as paper (Bunny, Spot, Pig, StanfordBunny, Armadillo, MaxPlanck, HappyBuddha)
- [ ] 1000 random vertex pairs per mesh
- [ ] Compare: path length, iteration count, wall time
- [ ] Verify: GFR path ≤ FlipOut path (within floating point tolerance)

**Output**: Results table matching paper format

---

## File Structure

```
geometry-central/
├── include/geometrycentral/surface/
│   └── funnel_geodesics.h
├── src/surface/
│   └── funnel_geodesics.cpp
└── test/src/
    └── funnel_geodesics_test.cpp
```

## Estimated LOC per Phase

| Phase | New LOC | Cumulative |
|-------|---------|------------|
| 0 | 20 | 20 |
| 1 | 80 | 100 |
| 2 | 40 | 140 |
| 3 | 120 | 260 |
| 4 | 60 | 320 |
| 5 | 50 | 370 |
| 6 | 100 | 470 |
| 7 | 60 | 530 |
| 8 | 120 | 650 |
| 9 | 50 | 700 |

**Total: ~700 LOC** (excluding tests)

## Risk Mitigation

1. **Numerical precision**: Use `double` throughout (not `float`)
2. **Boundary cases**: Test mesh boundaries explicitly in Phase 6
3. **Degenerate geometry**: Handle collinear vertices in funnel
4. **Performance**: Don't optimize until Phase 10 validates correctness
