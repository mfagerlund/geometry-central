#pragma once

// Funnel Geodesics (GFR) - Greedy Funnel Refinement
//
// An alternative to FlipOut for computing geodesic paths on triangle meshes.
// Uses sleeve refinement and quality-gated corner flipping instead of
// intrinsic edge flipping.
//
// Reference: See docs/GFR_PORT_PLAN.md for algorithm details.
//
// Ported from C# implementation in Colonel.Meshing.GreedyFunnelRefinement

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_point.h"

#include <memory>
#include <vector>

namespace geometrycentral {
namespace surface {

// Forward declarations
class FunnelGeodesicPath;

// ============================================================================
// Main API
// ============================================================================

/// Compute a geodesic path between two vertices using Funnel Refinement.
/// This is the main entry point - mirrors FlipEdgeNetwork::constructFromDijkstraPath()
std::unique_ptr<FunnelGeodesicPath> computeFunnelGeodesic(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    Vertex startVert,
    Vertex endVert);

/// Get cache statistics for the VeryDiscreteGeodesic pathfinder (for debugging)
struct CacheStats {
  size_t hits = 0;
  size_t misses = 0;
  size_t cacheSize = 0;
};
CacheStats getCacheStats();

/// Get timing breakdown for profiling
struct TimingStats {
  double aStarMs = 0;
  double flattenMs = 0;
  double straightenMs = 0;
  size_t pathCount = 0;
};
TimingStats getTimingStats();
void resetTimingStats();

// ============================================================================
// Result class
// ============================================================================

class FunnelGeodesicPath {
public:
  FunnelGeodesicPath(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom);

  // Query results
  double length() const;
  const std::vector<SurfacePoint>& getPath() const;

  // Statistics
  size_t iterationCount() const;
  size_t faceCount() const;  // faces in final sleeve

  // For visualization (optional)
  VertexPositionGeometry* posGeom = nullptr;
  std::vector<Vector3> getPathPolyline3D() const;

private:
  ManifoldSurfaceMesh& mesh;
  VertexPositionGeometry& geom;

  // Results
  std::vector<SurfacePoint> pathPoints;
  double pathLength = 0.0;
  size_t nIterations = 0;
  size_t nFaces = 0;

  // Internal: Sleeve data
  std::vector<Face> sleeveFaces;

  friend std::unique_ptr<FunnelGeodesicPath> computeFunnelGeodesic(
      ManifoldSurfaceMesh&, VertexPositionGeometry&, Vertex, Vertex);
};

// ============================================================================
// Internal structures (exposed for testing)
// ============================================================================

namespace funnel_internal {

/// 2D position of a vertex in the unfolded sleeve
struct FlatVertex {
  Vertex vertex;
  Vector2 pos;
};

/// Portal between two adjacent faces in the sleeve
struct Portal {
  Vector2 left, right;
  Vertex leftVert, rightVert;
};

/// Result of funnel algorithm
struct FunnelResult {
  std::vector<Vector2> waypoints2D;
  std::vector<Vertex> waypointVertices;
  double distance;
};

/// Corner that may want to flip
struct WaypointCorner {
  Vertex vertex;
  size_t faceIndex;
  double angleErrorDeg;
  // C#: public bool WantToFlip() => AngleErrorDeg > 0 && CanStraighten;
  // Note: CanStraighten is checked in computeFlipAction (canFlip)
  bool wantsToFlip() const { return angleErrorDeg > 0; }
};

/// Action to flip a corner (add/remove faces)
struct FlipAction {
  std::vector<Face> removeFaces;
  std::vector<Face> addFaces;
  size_t spliceAfterIndex;  // Insert addFaces after this face index
  size_t spliceBeforeIndex; // And before this face index
  bool canFlip;
};

// Phase 1: Flatten face strip to 2D
VertexData<Vector2> flattenSleeve(
    const std::vector<Face>& faces,
    Vertex entry,
    VertexPositionGeometry& geom);

// Phase 2: Build portals from flattened sleeve
std::vector<Portal> buildPortals(
    const std::vector<Face>& faces,
    const VertexData<Vector2>& flatPos);

// Phase 3: Run funnel algorithm
FunnelResult runFunnel(
    const std::vector<Portal>& portals,
    Vector2 entry,
    Vector2 exit);

// Phase 4: Build face strip from Dijkstra path
std::vector<Face> buildFaceStrip(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    Vertex start,
    Vertex end);

// Phase 4b: Build face strip via VeryDiscreteGeodesic (better initial corridors)
std::vector<Face> buildFaceStripVeryDiscrete(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    Vertex start,
    Vertex end);

// Phase 5: Analyze corners
std::vector<WaypointCorner> analyzeCorners(
    const std::vector<Face>& faces,
    const FunnelResult& funnel,
    const VertexData<Vector2>& flatPos);

// Phase 6: Compute flip action
FlipAction computeFlipAction(
    const std::vector<Face>& faces,
    const WaypointCorner& corner,
    ManifoldSurfaceMesh& mesh);

// Phase 7: Apply flip (returns new face list)
std::vector<Face> applyFlip(
    const std::vector<Face>& faces,
    const FlipAction& action);

} // namespace funnel_internal

} // namespace surface
} // namespace geometrycentral
