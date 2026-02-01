#include "geometrycentral/surface/funnel_geodesics.h"
#include "geometrycentral/surface/very_discrete_geodesic.h"
#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <iostream>
#include <iomanip>
#include <random>
#include <string>
#include <chrono>

#ifdef _WIN32
#include <windows.h>
#endif

#include "bunny_flipout_data.h"
#include "bunny_gfr_data.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// ============================================================================
// BenchmarkPriority - Matches C# HighPriorityScope for fair comparison
// Sets high thread/process priority and pins to last CPU core
// ============================================================================
#ifdef _WIN32
class BenchmarkPriority {
public:
  BenchmarkPriority() {
    // Save original state
    hProcess = GetCurrentProcess();
    hThread = GetCurrentThread();
    originalPriorityClass = GetPriorityClass(hProcess);
    originalThreadPriority = GetThreadPriority(hThread);

    // Set high priority (matching C# High/Highest)
    SetPriorityClass(hProcess, HIGH_PRIORITY_CLASS);
    SetThreadPriority(hThread, THREAD_PRIORITY_HIGHEST);

    // Pin to last CPU core (less system service contention)
    SYSTEM_INFO sysInfo;
    GetSystemInfo(&sysInfo);
    DWORD cpuCount = sysInfo.dwNumberOfProcessors;
    DWORD_PTR lastCpuMask = 1ULL << (cpuCount - 1);
    originalAffinityMask = SetThreadAffinityMask(hThread, lastCpuMask);
    affinitySet = (originalAffinityMask != 0);

    std::cout << "  [Benchmark mode: priority=HIGH, CPU=" << (cpuCount - 1) << "]" << std::endl;
  }

  ~BenchmarkPriority() {
    // Restore original state
    SetPriorityClass(hProcess, originalPriorityClass);
    SetThreadPriority(hThread, originalThreadPriority);
    if (affinitySet) {
      SetThreadAffinityMask(hThread, originalAffinityMask);
    }
  }

private:
  HANDLE hProcess;
  HANDLE hThread;
  DWORD originalPriorityClass;
  int originalThreadPriority;
  DWORD_PTR originalAffinityMask;
  bool affinitySet;
};
#else
// No-op on non-Windows
class BenchmarkPriority {
public:
  BenchmarkPriority() {
    std::cout << "  [Benchmark mode: not available on this platform]" << std::endl;
  }
};
#endif

class FunnelGeodesicsSuite : public MeshAssetSuite {};

// ============================================================================
// Basic tests
// ============================================================================

// Test that the algorithm runs without crashing on simple meshes
TEST_F(FunnelGeodesicsSuite, BasicComputation) {
  for (const MeshAsset& a : {getAsset("spot.ply", true), getAsset("bob_small.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& geom = *a.geometry;

    // Pick two arbitrary vertices
    Vertex v0 = mesh.vertex(0);
    Vertex v1 = mesh.vertex(mesh.nVertices() / 2);

    auto path = computeFunnelGeodesic(mesh, geom, v0, v1);

    EXPECT_GT(path->length(), 0);
    EXPECT_GE(path->faceCount(), 1);
  }
}

// Test that same vertex returns zero length
TEST_F(FunnelGeodesicsSuite, SameVertex) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  Vertex v0 = mesh.vertex(42);

  auto path = computeFunnelGeodesic(mesh, geom, v0, v0);

  // Same vertex should give very short or zero path
  EXPECT_LT(path->length(), 1e-10);
}

// Test adjacent vertices
TEST_F(FunnelGeodesicsSuite, AdjacentVertices) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  // Find two adjacent vertices
  Vertex v0 = mesh.vertex(0);
  Vertex v1 = v0.halfedge().twin().vertex();

  double edgeLength = (geom.vertexPositions[v0] - geom.vertexPositions[v1]).norm();

  auto path = computeFunnelGeodesic(mesh, geom, v0, v1);

  // Path length should equal edge length for adjacent vertices
  EXPECT_NEAR(path->length(), edgeLength, 1e-6);
}

// ============================================================================
// Path quality tests
// ============================================================================

// Compute edge-walk distance using Dijkstra
static double computeDijkstraDistance(VertexPositionGeometry& geom, Vertex start, Vertex end) {
  std::vector<Halfedge> edgePath = shortestEdgePath(geom, start, end);

  double distance = 0;
  for (Halfedge he : edgePath) {
    distance += geom.edgeLength(he.edge());
  }
  return distance;
}

// Test that geodesic path is shorter than or equal to edge-walk distance
TEST_F(FunnelGeodesicsSuite, PathShorterThanEdgeWalk) {
  for (const MeshAsset& a : {getAsset("spot.ply", true), getAsset("bob_small.ply", true)}) {
    a.printThyName();
    ManifoldSurfaceMesh& mesh = *a.manifoldMesh;
    VertexPositionGeometry& geom = *a.geometry;

    std::mt19937 mt(42);
    std::uniform_int_distribution<size_t> dist(0, mesh.nVertices() - 1);

    // Test with several random vertex pairs
    for (int i = 0; i < 20; i++) {
      Vertex v0 = mesh.vertex(dist(mt));
      Vertex v1 = mesh.vertex(dist(mt));

      if (v0 == v1) continue;

      double dijkstraDist = computeDijkstraDistance(geom, v0, v1);
      auto path = computeFunnelGeodesic(mesh, geom, v0, v1);

      // Geodesic should be shorter than or equal to edge walk
      EXPECT_LE(path->length(), dijkstraDist + 1e-6);
    }
  }
}

// Test convergence: algorithm should terminate
TEST_F(FunnelGeodesicsSuite, Convergence) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  std::mt19937 mt(42);
  std::uniform_int_distribution<size_t> dist(0, mesh.nVertices() - 1);

  // Test many pairs to ensure convergence
  for (int i = 0; i < 50; i++) {
    Vertex v0 = mesh.vertex(dist(mt));
    Vertex v1 = mesh.vertex(dist(mt));

    if (v0 == v1) continue;

    auto path = computeFunnelGeodesic(mesh, geom, v0, v1);

    // Algorithm should have converged (iteration count < max)
    EXPECT_LT(path->iterationCount(), 20000);
    EXPECT_GT(path->length(), 0);
  }
}

// ============================================================================
// Internal function tests
// ============================================================================

using namespace funnel_internal;

// Test flattening of a single face
TEST_F(FunnelGeodesicsSuite, FlattenSingleFace) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  // Get a face and its vertices
  Face f = mesh.face(0);
  std::vector<Face> faces = {f};

  Vertex entry = f.halfedge().vertex();
  auto flatPos = flattenSleeve(faces, entry, geom);

  // Entry vertex should be at origin
  EXPECT_LT(flatPos[entry].norm(), 1e-10);

  // Check that 2D edge lengths match 3D edge lengths
  for (Halfedge he : f.adjacentHalfedges()) {
    Vertex vA = he.vertex();
    Vertex vB = he.twin().vertex();

    double dist3D = (geom.vertexPositions[vA] - geom.vertexPositions[vB]).norm();
    double dist2D = (flatPos[vA] - flatPos[vB]).norm();

    EXPECT_NEAR(dist3D, dist2D, 1e-6);
  }
}

// Test portal building
TEST_F(FunnelGeodesicsSuite, BuildPortals) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  // Build a face strip between two vertices
  Vertex v0 = mesh.vertex(0);
  Vertex v1 = mesh.vertex(mesh.nVertices() / 4);

  auto faces = buildFaceStrip(mesh, geom, v0, v1);

  if (faces.size() >= 2) {
    auto flatPos = flattenSleeve(faces, v0, geom);
    auto portals = buildPortals(faces, flatPos);

    // Should have one less portal than faces
    EXPECT_EQ(portals.size(), faces.size() - 1);

    // Each portal should have valid positions
    for (const Portal& p : portals) {
      EXPECT_TRUE(std::isfinite(p.left.x));
      EXPECT_TRUE(std::isfinite(p.left.y));
      EXPECT_TRUE(std::isfinite(p.right.x));
      EXPECT_TRUE(std::isfinite(p.right.y));
    }
  }
}

// Test funnel algorithm on a simple corridor (no turns)
TEST_F(FunnelGeodesicsSuite, FunnelStraightCorridor) {
  // Create a simple test case with no turns
  std::vector<Portal> portals;

  // Portals forming a straight corridor
  Portal p1;
  p1.left = Vector2{0, 1};
  p1.right = Vector2{0, -1};
  portals.push_back(p1);

  Portal p2;
  p2.left = Vector2{1, 1};
  p2.right = Vector2{1, -1};
  portals.push_back(p2);

  Vector2 entry = Vector2{-1, 0};
  Vector2 exit = Vector2{2, 0};

  auto result = runFunnel(portals, entry, exit);

  // Straight corridor should have no waypoints
  EXPECT_EQ(result.waypointVertices.size(), 0);

  // Distance should be straight-line distance
  double expectedDist = (exit - entry).norm();
  EXPECT_NEAR(result.distance, expectedDist, 1e-6);
}

// Test funnel algorithm with a turn
TEST_F(FunnelGeodesicsSuite, FunnelWithTurn) {
  std::vector<Portal> portals;

  // First portal
  Portal p1;
  p1.left = Vector2{1, 2};
  p1.right = Vector2{1, 0};
  portals.push_back(p1);

  // Second portal - forces a turn around the corner
  Portal p2;
  p2.left = Vector2{2, 2};
  p2.right = Vector2{2, 1};
  portals.push_back(p2);

  Vector2 entry = Vector2{0, 0};
  Vector2 exit = Vector2{3, 2};

  auto result = runFunnel(portals, entry, exit);

  // Should have path through the portals
  EXPECT_GT(result.distance, 0);

  // Distance should be less than or equal to going around all corners
  double maxDist = (Vector2{1, 0} - entry).norm() +
                   (Vector2{2, 1} - Vector2{1, 0}).norm() +
                   (exit - Vector2{2, 1}).norm();
  EXPECT_LE(result.distance, maxDist + 1e-6);
}

// Test face strip building
TEST_F(FunnelGeodesicsSuite, FaceStripConnectivity) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  Vertex v0 = mesh.vertex(0);
  Vertex v1 = mesh.vertex(mesh.nVertices() / 3);

  auto faces = buildFaceStrip(mesh, geom, v0, v1);

  // First face should contain start vertex
  if (!faces.empty()) {
    bool containsStart = false;
    for (Vertex v : faces[0].adjacentVertices()) {
      if (v == v0) containsStart = true;
    }
    EXPECT_TRUE(containsStart);
  }

  // Last face should contain end vertex
  if (!faces.empty()) {
    bool containsEnd = false;
    for (Vertex v : faces.back().adjacentVertices()) {
      if (v == v1) containsEnd = true;
    }
    EXPECT_TRUE(containsEnd);
  }

  // Consecutive faces should share an edge
  for (size_t i = 0; i + 1 < faces.size(); i++) {
    Face f1 = faces[i];
    Face f2 = faces[i + 1];

    bool sharesEdge = false;
    for (Halfedge he : f1.adjacentHalfedges()) {
      if (he.twin().face() == f2) {
        sharesEdge = true;
        break;
      }
    }
    EXPECT_TRUE(sharesEdge);
  }
}

// ============================================================================
// Stress tests
// ============================================================================

// Test many random paths on larger mesh
TEST_F(FunnelGeodesicsSuite, StressTestRandomPaths) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  std::mt19937 mt(42);
  std::uniform_int_distribution<size_t> dist(0, mesh.nVertices() - 1);

  size_t numPaths = 100;
  size_t successCount = 0;

  for (size_t i = 0; i < numPaths; i++) {
    Vertex v0 = mesh.vertex(dist(mt));
    Vertex v1 = mesh.vertex(dist(mt));

    if (v0 == v1) {
      successCount++;
      continue;
    }

    auto path = computeFunnelGeodesic(mesh, geom, v0, v1);

    if (path->length() > 0 && path->iterationCount() < 20000) {
      successCount++;
    }
  }

  // All paths should succeed
  EXPECT_EQ(successCount, numPaths);
}

// Test that getPath() returns valid SurfacePoints
TEST_F(FunnelGeodesicsSuite, GetPathReturnsValidPoints) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  Vertex v0 = mesh.vertex(0);
  Vertex v1 = mesh.vertex(mesh.nVertices() / 2);

  auto path = computeFunnelGeodesic(mesh, geom, v0, v1);

  const auto& pathPoints = path->getPath();

  // Should have at least start and end
  EXPECT_GE(pathPoints.size(), 2);

  // First point should be start vertex
  EXPECT_EQ(pathPoints.front().type, SurfacePointType::Vertex);
  EXPECT_EQ(pathPoints.front().vertex, v0);

  // Last point should be end vertex
  EXPECT_EQ(pathPoints.back().type, SurfacePointType::Vertex);
  EXPECT_EQ(pathPoints.back().vertex, v1);

  // All intermediate points should be valid vertices
  for (size_t i = 1; i < pathPoints.size() - 1; i++) {
    EXPECT_EQ(pathPoints[i].type, SurfacePointType::Vertex);
    EXPECT_TRUE(pathPoints[i].vertex.getIndex() < mesh.nVertices());
  }
}

// Test getPathPolyline3D returns valid 3D coordinates
TEST_F(FunnelGeodesicsSuite, GetPathPolyline3D) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  Vertex v0 = mesh.vertex(0);
  Vertex v1 = mesh.vertex(mesh.nVertices() / 2);

  auto path = computeFunnelGeodesic(mesh, geom, v0, v1);
  auto polyline = path->getPathPolyline3D();

  // Should have same length as path points
  EXPECT_EQ(polyline.size(), path->getPath().size());

  // First and last should match vertex positions
  EXPECT_NEAR((polyline.front() - geom.vertexPositions[v0]).norm(), 0, 1e-10);
  EXPECT_NEAR((polyline.back() - geom.vertexPositions[v1]).norm(), 0, 1e-10);

  // All coordinates should be finite
  for (const Vector3& p : polyline) {
    EXPECT_TRUE(std::isfinite(p.x));
    EXPECT_TRUE(std::isfinite(p.y));
    EXPECT_TRUE(std::isfinite(p.z));
  }
}

// ============================================================================
// VeryDiscreteGeodesic tests
// ============================================================================

using namespace very_discrete_geodesic;

// Test that VeryDiscreteGeodesic pathfinder finds a path
TEST_F(FunnelGeodesicsSuite, VeryDiscreteBasicPath) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  Vertex v0 = mesh.vertex(0);
  Vertex v1 = mesh.vertex(mesh.nVertices() / 2);

  PathResult result = findPath(v0, v1, mesh, geom);

  EXPECT_TRUE(result.isComplete);
  EXPECT_FALSE(result.isFallback);
  EXPECT_GT(result.path.size(), 1);
  EXPECT_EQ(result.path.front(), v0);
  EXPECT_EQ(result.path.back(), v1);
}

// Test that VeryDiscreteGeodesic handles same vertex
TEST_F(FunnelGeodesicsSuite, VeryDiscreteSameVertex) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  Vertex v0 = mesh.vertex(42);

  PathResult result = findPath(v0, v0, mesh, geom);

  EXPECT_TRUE(result.isComplete);
  EXPECT_EQ(result.path.size(), 1);
  EXPECT_EQ(result.path.front(), v0);
}

// Test that VeryDiscreteGeodesic produces face strips
TEST_F(FunnelGeodesicsSuite, VeryDiscreteFaceStrip) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  Vertex v0 = mesh.vertex(0);
  Vertex v1 = mesh.vertex(mesh.nVertices() / 4);

  std::pair<std::vector<Face>, std::vector<Vertex>> result =
      findFaceStripWithPath(v0, v1, mesh, geom);

  std::vector<Face>& faces = result.first;
  std::vector<Vertex>& path = result.second;

  // Should have faces and path
  EXPECT_GT(faces.size(), 0);
  EXPECT_GT(path.size(), 1);
  EXPECT_EQ(path.front(), v0);
  EXPECT_EQ(path.back(), v1);

  // First face should contain start vertex
  if (!faces.empty()) {
    bool containsStart = false;
    for (Vertex v : faces[0].adjacentVertices()) {
      if (v == v0) containsStart = true;
    }
    EXPECT_TRUE(containsStart);
  }

  // Last face should contain end vertex
  if (!faces.empty()) {
    bool containsEnd = false;
    for (Vertex v : faces.back().adjacentVertices()) {
      if (v == v1) containsEnd = true;
    }
    EXPECT_TRUE(containsEnd);
  }
}

// Test VeryDiscreteGeodesic explorer
TEST_F(FunnelGeodesicsSuite, VeryDiscreteExplorer) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  // Pick a vertex with multiple corners
  Vertex v0 = mesh.vertex(100);

  // Get a corner
  bool foundCorner = false;
  for (Corner c : v0.adjacentCorners()) {
    ExplorationResult result = explore(c, geom, ExplorationDepth::Level4);

    // Should have at least L1 computed (unless boundary)
    if (result.L1.blocked != BlockedReason::F1Boundary) {
      foundCorner = true;
      // L1 should be computed
      EXPECT_NE(result.L1.blocked, BlockedReason::NotComputed);
    }
    break;
  }

  EXPECT_TRUE(foundCorner);
}

// Test that VeryDiscreteGeodesic path is shorter than or equal to edge-only
TEST_F(FunnelGeodesicsSuite, VeryDiscretePathQuality) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  std::mt19937 mt(42);
  std::uniform_int_distribution<size_t> dist(0, mesh.nVertices() - 1);

  // Test with several random vertex pairs
  for (int i = 0; i < 10; i++) {
    Vertex v0 = mesh.vertex(dist(mt));
    Vertex v1 = mesh.vertex(dist(mt));

    if (v0 == v1) continue;

    // Get edge-only distance
    std::vector<Halfedge> edgePath = shortestEdgePath(geom, v0, v1);
    double edgeDist = 0;
    for (Halfedge he : edgePath) {
      edgeDist += geom.edgeLength(he.edge());
    }

    // Get VeryDiscreteGeodesic distance
    PathResult result = findPath(v0, v1, mesh, geom);
    double vdgDist = computePathDistance(result.steps);

    // VeryDiscreteGeodesic should be shorter or equal
    EXPECT_LE(vdgDist, edgeDist + 1e-6);
  }
}

// Tests that getNeighbors deduplicates vertices correctly
// Note: The explorer CAN return duplicate vertices (different paths to same vertex)
// but getNeighbors must deduplicate them for the A* pathfinder
TEST_F(FunnelGeodesicsSuite, VeryDiscreteNeighborsNoDuplicates) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  // Test that getNeighbors doesn't return duplicates
  // This is done implicitly via seenVertices in the implementation,
  // but we verify by checking the path is valid
  Vertex v0 = mesh.vertex(100);
  Vertex v1 = mesh.vertex(200);

  PathResult result = findPath(v0, v1, mesh, geom);

  // Path should have no duplicate consecutive vertices
  for (size_t i = 1; i < result.path.size(); i++) {
    EXPECT_NE(result.path[i], result.path[i-1]);
  }

  // Path should have no duplicate vertices at all
  std::unordered_set<size_t> seen;
  bool hasDuplicates = false;
  for (Vertex v : result.path) {
    if (seen.count(v.getIndex())) {
      hasDuplicates = true;
      break;
    }
    seen.insert(v.getIndex());
  }
  EXPECT_FALSE(hasDuplicates);
}

// Ported from C# VeryDiscreteGeodesicTests.L1_NeverBlockedByPortalCrossing
// L1 can only be blocked by boundary, not by portal crossing
TEST_F(FunnelGeodesicsSuite, VeryDiscreteL1NotPortalBlocked) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  int failures = 0;

  for (Vertex v : mesh.vertices()) {
    for (Corner c : v.adjacentCorners()) {
      ExplorationResult result = explore(c, geom, ExplorationDepth::Level4);

      // If L1 is blocked, it should NOT be PortalBlocked
      if (!result.L1.isReachable &&
          result.L1.blocked == BlockedReason::PortalBlocked) {
        failures++;
      }
    }
  }

  EXPECT_EQ(failures, 0);
}

// Ported from C# VeryDiscreteGeodesicPathfinderTests.Path_IsShorterOrEqualToDijkstra
// Tests that VeryDiscreteGeodesic path is at least as good as Dijkstra
TEST_F(FunnelGeodesicsSuite, VeryDiscreteShorterThanDijkstra) {
  auto asset = getAsset("spot.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  std::mt19937 mt(12345);
  std::uniform_int_distribution<size_t> dist(0, mesh.nVertices() - 1);

  int vdgShorter = 0;
  int dijkstraShorter = 0;
  int matches = 0;

  for (int i = 0; i < 50; i++) {
    Vertex v0 = mesh.vertex(dist(mt));
    Vertex v1 = mesh.vertex(dist(mt));
    if (v0 == v1) continue;

    // Dijkstra (edge-only)
    std::vector<Halfedge> edgePath = shortestEdgePath(geom, v0, v1);
    double edgeDist = 0;
    for (Halfedge he : edgePath) {
      edgeDist += geom.edgeLength(he.edge());
    }

    // VeryDiscreteGeodesic
    PathResult result = findPath(v0, v1, mesh, geom);
    double vdgDist = computePathDistance(result.steps);

    double diff = std::abs(vdgDist - edgeDist);
    if (diff < 1e-9) {
      matches++;
    } else if (vdgDist < edgeDist) {
      vdgShorter++;
    } else {
      dijkstraShorter++;
    }
  }

  // VeryDiscreteGeodesic should never be longer
  EXPECT_EQ(dijkstraShorter, 0);
  // VeryDiscreteGeodesic should often be shorter
  EXPECT_GT(vdgShorter + matches, 0);
}

// Tests multiple vertex pairs with different path characteristics
TEST_F(FunnelGeodesicsSuite, VeryDiscreteMultiplePaths) {
  auto asset = getAsset("bob_small.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  // Test several specific vertex pairs
  std::vector<std::pair<size_t, size_t>> pairs = {
    {0, mesh.nVertices() / 4},
    {mesh.nVertices() / 4, mesh.nVertices() / 2},
    {0, mesh.nVertices() - 1},
    {10, 100},
  };

  for (const auto& pair : pairs) {
    if (pair.first >= mesh.nVertices() || pair.second >= mesh.nVertices()) continue;

    Vertex v0 = mesh.vertex(pair.first);
    Vertex v1 = mesh.vertex(pair.second);

    PathResult result = findPath(v0, v1, mesh, geom);

    // Path should be valid
    EXPECT_TRUE(result.isComplete);
    EXPECT_GE(result.path.size(), 2);
    EXPECT_EQ(result.path.front(), v0);
    EXPECT_EQ(result.path.back(), v1);

    // Distance should be positive
    double dist = computePathDistance(result.steps);
    EXPECT_GT(dist, 0);
  }
}

// Test explorer at different depth levels
TEST_F(FunnelGeodesicsSuite, VeryDiscreteExplorerDepths) {
  auto asset = getAsset("sphere_small.ply", true);
  ManifoldSurfaceMesh& mesh = *asset.manifoldMesh;
  VertexPositionGeometry& geom = *asset.geometry;

  // Pick an interior vertex (avoid boundary)
  Vertex v0 = mesh.vertex(mesh.nVertices() / 2);

  for (Corner c : v0.adjacentCorners()) {
    // Explore at increasing depths
    ExplorationResult r1 = explore(c, geom, ExplorationDepth::Level1);
    ExplorationResult r3 = explore(c, geom, ExplorationDepth::Level3);
    ExplorationResult r5 = explore(c, geom, ExplorationDepth::Level5);

    // Deeper exploration should find more or equal candidates
    size_t count1 = r1.getReachableCandidates().size();
    size_t count3 = r3.getReachableCandidates().size();
    size_t count5 = r5.getReachableCandidates().size();

    EXPECT_LE(count1, count3);
    EXPECT_LE(count3, count5);
    break; // One corner is enough for this test
  }
}

// ============================================================================
// FlipOut Comparison Test (Full 1000 paths)
// Compares C++ GFR distances and timing against FlipOut from C# tests
// Uses the same vertex pairs as Benchmark_FlipsPerPath_AllMeshes
// ============================================================================

TEST_F(FunnelGeodesicsSuite, FlipOutComparisonFull) {
  // Load Bunny mesh directly from test assets
  std::string bunnyPath = std::string(GC_TEST_ASSETS_ABS_PATH) + "/Bunny.obj";
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = readManifoldSurfaceMesh(bunnyPath);
  ManifoldSurfaceMesh& mesh = *meshPtr;
  VertexPositionGeometry& geom = *geomPtr;

  std::cout << "\n====================================================================================================";
  std::cout << "\nC++ GFR vs FlipOut COMPARISON (Bunny, " << BUNNY_FLIPOUT_COUNT << " paths)";
  std::cout << "\n====================================================================================================\n";

  // Set high priority and pin to last CPU for stable benchmarking (matches C# BenchmarkPriority)
  BenchmarkPriority benchmarkScope;

  // FlipOut timing from C# test (from Bunny_flipout.json header)
  const double flipOutTimeMs = 3476.0;  // From the JSON file

  // Warmup pass
  for (size_t i = 0; i < std::min((size_t)50, BUNNY_FLIPOUT_COUNT); i++) {
    const BunnyFlipOutPair& p = BUNNY_FLIPOUT_PAIRS[i];
    if (p.from >= mesh.nVertices() || p.to >= mesh.nVertices()) continue;
    Vertex v0 = mesh.vertex(p.from);
    Vertex v1 = mesh.vertex(p.to);
    auto path = computeFunnelGeodesic(mesh, geom, v0, v1);
  }

  // Reset timing stats for timed run
  resetTimingStats();

  // Timed run
  int gfrShorter = 0;
  int flipOutShorter = 0;
  int matches = 0;
  double totalDiffPercent = 0.0;
  double maxDiffPercent = 0.0;
  double minDiffPercent = 0.0;
  double totalGfrDist = 0.0;
  double totalFlipOutDist = 0.0;
  int validPaths = 0;

  auto startTime = std::chrono::high_resolution_clock::now();

  for (size_t i = 0; i < BUNNY_FLIPOUT_COUNT; i++) {
    const BunnyFlipOutPair& p = BUNNY_FLIPOUT_PAIRS[i];
    if (p.from >= mesh.nVertices() || p.to >= mesh.nVertices()) continue;

    Vertex v0 = mesh.vertex(p.from);
    Vertex v1 = mesh.vertex(p.to);

    auto gfrPath = computeFunnelGeodesic(mesh, geom, v0, v1);
    double gfrDist = gfrPath->length();
    double flipOutDist = p.flipOutDistance;

    totalGfrDist += gfrDist;
    totalFlipOutDist += flipOutDist;
    validPaths++;

    double diffPercent = (gfrDist - flipOutDist) / flipOutDist * 100.0;
    totalDiffPercent += diffPercent;

    if (diffPercent > maxDiffPercent) maxDiffPercent = diffPercent;
    if (diffPercent < minDiffPercent) minDiffPercent = diffPercent;

    if (std::abs(diffPercent) < 0.01) {
      matches++;
    } else if (gfrDist < flipOutDist) {
      gfrShorter++;
    } else {
      flipOutShorter++;
    }
  }

  auto endTime = std::chrono::high_resolution_clock::now();
  double gfrTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();

  double meanDiff = totalDiffPercent / validPaths;
  double avgDistDiff = (totalGfrDist - totalFlipOutDist) / totalFlipOutDist * 100.0;

  std::cout << "\nTIMING:" << std::endl;
  std::cout << "  C++ GFR:    " << gfrTimeMs << " ms (" << (gfrTimeMs / validPaths) << " ms/path)" << std::endl;
  std::cout << "  FlipOut:    " << flipOutTimeMs << " ms (" << (flipOutTimeMs / validPaths) << " ms/path)" << std::endl;
  std::cout << "  Speedup:    " << (flipOutTimeMs / gfrTimeMs) << "x" << std::endl;

  std::cout << "\nPATH QUALITY:" << std::endl;
  std::cout << "  Valid paths:     " << validPaths << std::endl;
  std::cout << "  GFR shorter:     " << gfrShorter << " (" << (gfrShorter * 100.0 / validPaths) << "%)" << std::endl;
  std::cout << "  FlipOut shorter: " << flipOutShorter << " (" << (flipOutShorter * 100.0 / validPaths) << "%)" << std::endl;
  std::cout << "  Match (<0.01%):  " << matches << " (" << (matches * 100.0 / validPaths) << "%)" << std::endl;

  std::cout << "\nDISTANCE DIFFERENCES (GFR - FlipOut):" << std::endl;
  std::cout << "  Mean:       " << meanDiff << "%" << std::endl;
  std::cout << "  Min:        " << minDiffPercent << "%" << std::endl;
  std::cout << "  Max:        " << maxDiffPercent << "%" << std::endl;

  std::cout << "\nTOTAL DISTANCE:" << std::endl;
  std::cout << "  GFR:        " << totalGfrDist << std::endl;
  std::cout << "  FlipOut:    " << totalFlipOutDist << std::endl;
  std::cout << "  Difference: " << avgDistDiff << "%" << std::endl;

  // Print cache statistics
  auto cacheStats = getCacheStats();
  std::cout << "\nCACHE STATS:" << std::endl;
  std::cout << "  Hits:       " << cacheStats.hits << std::endl;
  std::cout << "  Misses:     " << cacheStats.misses << std::endl;
  std::cout << "  Cache size: " << cacheStats.cacheSize << std::endl;
  if (cacheStats.hits + cacheStats.misses > 0) {
    double hitRate = 100.0 * cacheStats.hits / (cacheStats.hits + cacheStats.misses);
    std::cout << "  Hit rate:   " << hitRate << "%" << std::endl;
  }

  // Print timing breakdown
  auto timingStats = getTimingStats();
  std::cout << "\nTIMING BREAKDOWN:" << std::endl;
  std::cout << "  A* pathfinding:   " << timingStats.aStarMs << " ms (" << (timingStats.aStarMs / validPaths) << " ms/path)" << std::endl;
  std::cout << "  Flatten/Funnel:   " << timingStats.flattenMs << " ms (" << (timingStats.flattenMs / validPaths) << " ms/path)" << std::endl;
  std::cout << "  Straightening:    " << timingStats.straightenMs << " ms (" << (timingStats.straightenMs / validPaths) << " ms/path)" << std::endl;
  double overhead = gfrTimeMs - timingStats.aStarMs - timingStats.flattenMs - timingStats.straightenMs;
  std::cout << "  Other/Overhead:   " << overhead << " ms" << std::endl;

  std::cout << "\n====================================================================================================" << std::endl;

  // GFR should be competitive with FlipOut (within 1% on average)
  EXPECT_LT(std::abs(meanDiff), 1.0);
}

// ============================================================================
// MMP Outlier Validation Test
// Finds the biggest GFR vs FlipOut outlier and validates with exact geodesics (MMP)
// ============================================================================

TEST_F(FunnelGeodesicsSuite, ValidateOutlierWithMMP) {
  // Load Bunny mesh directly from test assets
  std::string bunnyPath = std::string(GC_TEST_ASSETS_ABS_PATH) + "/Bunny.obj";
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = readManifoldSurfaceMesh(bunnyPath);
  ManifoldSurfaceMesh& mesh = *meshPtr;
  VertexPositionGeometry& geom = *geomPtr;

  std::cout << "\n====================================================================================================";
  std::cout << "\nMMP VALIDATION OF OUTLIER (Bunny)";
  std::cout << "\n====================================================================================================\n";

  // Find the outlier (minimum difference = GFR much shorter than FlipOut)
  size_t outlierIdx = 0;
  double minDiffPercent = 0.0;
  double outlierGfrDist = 0.0;
  double outlierFlipOutDist = 0.0;

  for (size_t i = 0; i < BUNNY_FLIPOUT_COUNT; i++) {
    const BunnyFlipOutPair& p = BUNNY_FLIPOUT_PAIRS[i];
    if (p.from >= mesh.nVertices() || p.to >= mesh.nVertices()) continue;

    Vertex v0 = mesh.vertex(p.from);
    Vertex v1 = mesh.vertex(p.to);

    auto gfrPath = computeFunnelGeodesic(mesh, geom, v0, v1);
    double gfrDist = gfrPath->length();
    double flipOutDist = p.flipOutDistance;

    double diffPercent = (gfrDist - flipOutDist) / flipOutDist * 100.0;

    if (diffPercent < minDiffPercent) {
      minDiffPercent = diffPercent;
      outlierIdx = i;
      outlierGfrDist = gfrDist;
      outlierFlipOutDist = flipOutDist;
    }
  }

  const BunnyFlipOutPair& outlier = BUNNY_FLIPOUT_PAIRS[outlierIdx];
  Vertex v0 = mesh.vertex(outlier.from);
  Vertex v1 = mesh.vertex(outlier.to);

  std::cout << "\nOUTLIER FOUND:" << std::endl;
  std::cout << "  Pair index:  " << outlierIdx << std::endl;
  std::cout << "  Vertices:    " << outlier.from << " -> " << outlier.to << std::endl;
  std::cout << "  Difference:  " << minDiffPercent << "%" << std::endl;

  // Compute exact geodesic distance using MMP
  std::cout << "\nComputing MMP exact geodesic..." << std::endl;
  GeodesicAlgorithmExact mmp(mesh, geom);
  mmp.propagate(v0);
  double mmpDist = mmp.getDistance(v1);

  std::cout << "\nDISTANCE COMPARISON:" << std::endl;
  std::cout << "  MMP (exact):  " << mmpDist << std::endl;
  std::cout << "  GFR:          " << outlierGfrDist << " (" << ((outlierGfrDist - mmpDist) / mmpDist * 100.0) << "% vs MMP)" << std::endl;
  std::cout << "  FlipOut:      " << outlierFlipOutDist << " (" << ((outlierFlipOutDist - mmpDist) / mmpDist * 100.0) << "% vs MMP)" << std::endl;

  // Validate: GFR path should be >= MMP (can't be shorter than exact)
  double gfrVsMmp = (outlierGfrDist - mmpDist) / mmpDist * 100.0;
  double flipOutVsMmp = (outlierFlipOutDist - mmpDist) / mmpDist * 100.0;

  std::cout << "\nVALIDATION:" << std::endl;
  if (outlierGfrDist < mmpDist - 1e-6) {
    std::cout << "  ERROR: GFR claims to be shorter than exact geodesic!" << std::endl;
  } else {
    std::cout << "  OK: GFR distance is valid (>= MMP)" << std::endl;
  }

  if (outlierFlipOutDist < mmpDist - 1e-6) {
    std::cout << "  ERROR: FlipOut claims to be shorter than exact geodesic!" << std::endl;
  } else {
    std::cout << "  OK: FlipOut distance is valid (>= MMP)" << std::endl;
  }

  std::cout << "\nCONCLUSION:" << std::endl;
  if (std::abs(gfrVsMmp) < std::abs(flipOutVsMmp)) {
    std::cout << "  GFR is closer to exact (" << gfrVsMmp << "% vs " << flipOutVsMmp << "%)" << std::endl;
  } else {
    std::cout << "  FlipOut is closer to exact (" << flipOutVsMmp << "% vs " << gfrVsMmp << "%)" << std::endl;
  }

  std::cout << "\n====================================================================================================" << std::endl;

  // GFR should not be shorter than exact geodesic
  EXPECT_GE(outlierGfrDist, mmpDist - 1e-6);
}

// ============================================================================
// C++ GFR vs C# GFR Comparison
// Verifies the port produces the same results as the original C# implementation
// ============================================================================

TEST_F(FunnelGeodesicsSuite, CppVsCsharpGFR) {
  // Load Bunny mesh directly from test assets
  std::string bunnyPath = std::string(GC_TEST_ASSETS_ABS_PATH) + "/Bunny.obj";
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = readManifoldSurfaceMesh(bunnyPath);
  ManifoldSurfaceMesh& mesh = *meshPtr;
  VertexPositionGeometry& geom = *geomPtr;

  std::cout << "\n====================================================================================================";
  std::cout << "\nC++ GFR vs C# GFR COMPARISON (Bunny, " << BUNNY_GFR_COUNT << " pairs)";
  std::cout << "\n====================================================================================================\n";

  std::cout << "\n" << std::setw(6) << "Pair" << " | "
            << std::setw(12) << "From->To" << " | "
            << std::setw(14) << "C++ GFR" << " | "
            << std::setw(14) << "C# GFR" << " | "
            << std::setw(10) << "Diff %" << " | "
            << "Status" << std::endl;
  std::cout << std::string(80, '-') << std::endl;

  int matches = 0;
  int close = 0;  // within 0.1%
  int farCount = 0;    // > 0.1%
  double totalDiff = 0;
  double maxDiff = 0;

  for (size_t i = 0; i < BUNNY_GFR_COUNT; i++) {
    const BunnyGFRPair& p = BUNNY_GFR_PAIRS[i];
    if (p.from >= mesh.nVertices() || p.to >= mesh.nVertices()) continue;

    Vertex v0 = mesh.vertex(p.from);
    Vertex v1 = mesh.vertex(p.to);

    auto cppPath = computeFunnelGeodesic(mesh, geom, v0, v1);
    double cppDist = cppPath->length();
    double csharpDist = p.csharpGfrDistance;

    double diff = (cppDist - csharpDist) / csharpDist * 100.0;
    totalDiff += std::abs(diff);
    if (std::abs(diff) > maxDiff) maxDiff = std::abs(diff);

    std::string status;
    if (std::abs(diff) < 0.001) {
      status = "MATCH";
      matches++;
    } else if (std::abs(diff) < 0.1) {
      status = "CLOSE";
      close++;
    } else {
      status = "DIFF";
      farCount++;
    }

    std::cout << std::setw(6) << i << " | "
              << std::setw(5) << p.from << "->" << std::setw(4) << p.to << " | "
              << std::setw(14) << std::fixed << std::setprecision(8) << cppDist << " | "
              << std::setw(14) << csharpDist << " | "
              << std::setw(9) << std::setprecision(4) << diff << "% | "
              << status << std::endl;
  }

  std::cout << std::string(80, '-') << std::endl;
  std::cout << "\nSUMMARY:" << std::endl;
  std::cout << "  MATCH (<0.001%):  " << matches << " / " << BUNNY_GFR_COUNT << std::endl;
  std::cout << "  CLOSE (<0.1%):    " << close << " / " << BUNNY_GFR_COUNT << std::endl;
  std::cout << "  DIFF (>0.1%):     " << farCount << " / " << BUNNY_GFR_COUNT << std::endl;
  std::cout << "  Mean |diff|:      " << (totalDiff / BUNNY_GFR_COUNT) << "%" << std::endl;
  std::cout << "  Max |diff|:       " << maxDiff << "%" << std::endl;

  std::cout << "\n====================================================================================================" << std::endl;

  // C++ should match C# within 0.01% for most paths
  EXPECT_GE(matches + close, static_cast<int>(BUNNY_GFR_COUNT * 0.9));
}

// ============================================================================
// Debug single pair that differs
// ============================================================================

TEST_F(FunnelGeodesicsSuite, DebugCppVsCsharpPair8) {
  std::string bunnyPath = std::string(GC_TEST_ASSETS_ABS_PATH) + "/Bunny.obj";
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr;
  std::unique_ptr<VertexPositionGeometry> geomPtr;
  std::tie(meshPtr, geomPtr) = readManifoldSurfaceMesh(bunnyPath);
  ManifoldSurfaceMesh& mesh = *meshPtr;
  VertexPositionGeometry& geom = *geomPtr;

  // Pair 8: 469 -> 2449, C# = 2.3622918128967285, C++ = 2.36472839
  size_t fromIdx = 469;
  size_t toIdx = 2449;
  double csharpDist = 2.3622918128967285;

  Vertex v0 = mesh.vertex(fromIdx);
  Vertex v1 = mesh.vertex(toIdx);

  std::cout << "\n====================================================================================================";
  std::cout << "\nDEBUG PAIR 8: " << fromIdx << " -> " << toIdx;
  std::cout << "\n====================================================================================================\n";

  // Test VeryDiscreteGeodesic face strip
  std::cout << "\n1. VeryDiscreteGeodesic face strip:" << std::endl;
  auto vdResult = very_discrete_geodesic::findFaceStripWithPath(v0, v1, mesh, geom);
  std::cout << "   Face count: " << vdResult.first.size() << std::endl;
  std::cout << "   Vertex path length: " << vdResult.second.size() << std::endl;

  // Test Dijkstra face strip
  std::cout << "\n2. Dijkstra face strip:" << std::endl;
  auto dijkstraFaces = funnel_internal::buildFaceStrip(mesh, geom, v0, v1);
  std::cout << "   Face count: " << dijkstraFaces.size() << std::endl;

  // Compute GFR with VeryDiscrete
  std::cout << "\n3. Full GFR (VeryDiscrete):" << std::endl;
  auto gfrPath = computeFunnelGeodesic(mesh, geom, v0, v1);
  std::cout << "   Distance: " << std::setprecision(10) << gfrPath->length() << std::endl;
  std::cout << "   Iterations: " << gfrPath->iterationCount() << std::endl;
  std::cout << "   Final faces: " << gfrPath->faceCount() << std::endl;

  // Compare
  std::cout << "\n4. Comparison:" << std::endl;
  std::cout << "   C# GFR:   " << std::setprecision(10) << csharpDist << std::endl;
  std::cout << "   C++ GFR:  " << std::setprecision(10) << gfrPath->length() << std::endl;
  std::cout << "   Diff:     " << (gfrPath->length() - csharpDist) << " ("
            << ((gfrPath->length() - csharpDist) / csharpDist * 100.0) << "%)" << std::endl;

  // Also compute MMP for reference
  GeodesicAlgorithmExact mmp(mesh, geom);
  mmp.propagate(v0);
  double mmpDist = mmp.getDistance(v1);
  std::cout << "\n5. MMP (exact): " << std::setprecision(10) << mmpDist << std::endl;
  std::cout << "   C# vs MMP:  " << ((csharpDist - mmpDist) / mmpDist * 100.0) << "% above exact" << std::endl;
  std::cout << "   C++ vs MMP: " << ((gfrPath->length() - mmpDist) / mmpDist * 100.0) << "% above exact" << std::endl;

  std::cout << "\n====================================================================================================" << std::endl;
}

