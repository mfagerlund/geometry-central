#include "geometrycentral/surface/funnel_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "load_test_meshes.h"

#include "gtest/gtest.h"

#include <iostream>
#include <random>
#include <string>

using namespace geometrycentral;
using namespace geometrycentral::surface;

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
